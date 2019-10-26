#include "graph.h"
#include <algorithm>
#include <mutex>
#include <cstring>
#include "tbb/scalable_allocator.h"
#include <stdio.h>
#include <stdlib.h>


std::mutex gpu_lock;

const char A_NT = 0;
const char C_NT = 1;
const char G_NT = 2;
const char T_NT = 3;
const char N_NT = 4;

int mat_offset[] = { 0, 1, 3, 6 };                                                             
int *d_sub_mat;

int err;                            
int check_status = 0;

char* d_ref_seq;
char* d_query_seq;
uint64_t* d_seed_offsets;
uint32_t* d_index_table;
uint64_t* d_pos_table;

__global__
void compress_string (uint32_t n, char* src_seq, char* dst_seq){ 
    int thread_id = threadIdx.x;
    int block_dim = blockDim.x;
    int grid_dim = blockDim.x;
    int block_id = blockIdx.x;

    int stride = block_dim * grid_dim;
    uint32_t start = block_dim * block_id + thread_id;

    uint32_t len = n;

    for (uint32_t i = start; i < len; i += stride) {
        char ch = src_seq[i];
        char dst = N_NT;
        if ((ch == 'a') || (ch == 'A'))
            dst = A_NT;
        else if ((ch == 'c') || (ch == 'C'))
            dst = C_NT;
        else if ((ch == 'g') || (ch == 'G'))
            dst = G_NT;
        else if ((ch == 't') || (ch == 'T'))
            dst = T_NT;
        dst_seq[i] = dst;
    }
}

__global__
void find_hit_offsets (int num_seeds, int query_start, char* d_ref_seq, char* d_query_seq, uint32_t* d_index_table, uint64_t* d_pos_table, uint64_t *seed_offsets, int *d_sub_mat){
    int thread_id = threadIdx.x;
    int block_dim = blockDim.x;
    int grid_dim = gridDim.x;
    int block_id = blockIdx.x;

    int id_start = block_dim * block_id;
    int stride = block_dim * grid_dim;

    uint32_t start, end;
    uint32_t q_start;
    uint64_t seed_offset;
    uint32_t seed;

    __shared__ uint32_t r_starts[4096];
    __shared__ uint32_t q_starts[4096];
    __shared__ uint32_t total[32];
    __shared__ uint32_t total_hits;
    __shared__ int y_drop;

    y_drop = 9430;

    for (int id = id_start; id < num_seeds; id += stride) {
        total[thread_id] = 0;

        // id+thread_id is a valid address
        if (id + thread_id < num_seeds) {
            seed_offset = seed_offsets[id+thread_id];

            seed = (seed_offset >> 32);
            q_start = ((seed_offset << 32) >> 32);

            // start and end from the seed block_id table
            end = d_index_table[seed];
            start = 0;
            if (seed > 0) {
                start = d_index_table[seed-1];
            }
        }
        else {
            start = 0;
            end = 0;
        }
        total[thread_id] += end-start;
        __syncthreads();
        
        // cumulative sum
        for(int k = 0; k < 5; k++){
            if(thread_id >= (1<<k)){
                total[thread_id] = total[thread_id] + total[thread_id-(1<<k)];
            }
            __syncthreads();
        }

        total_hits = total[31];
        int addr_start = (thread_id == 0) ? 0 : total[thread_id-1];
        for (uint32_t p = start; p < end; p++) { 
            if (p+addr_start-start < 4096) { 
                int index_el = p+addr_start-start;
                r_starts[index_el] = d_pos_table[p];
                q_starts[index_el] = q_start + query_start;
            }
        }
        __syncthreads();
    }

    int ref_loc;
    int query_loc;
    int seed_score;
    int right_ext_score;
    int left_ext_score;
    int ext_score;
    int match_score;
    int max_score;
    int total_score;

    for (int id = id_start; id < total_hits; id += stride) {
        ref_loc = r_starts[id+thread_id];
        query_loc = q_starts[id+thread_id];
        
        match_score = d_sub_mat[d_ref_seq[ref_loc]*5+d_query_seq[query_loc]];
        seed_score = match_score;
        
        for(int i = 1; i < 19; i++){
            match_score = d_sub_mat[d_ref_seq[ref_loc+i]*5+d_query_seq[query_loc+i]];
            seed_score = seed_score + match_score;
        }

        max_score = 0;
        ref_loc = r_starts[id+thread_id] + 19;
        query_loc = q_starts[id+thread_id] + 19;
        ext_score = d_sub_mat[d_ref_seq[ref_loc]*5+d_query_seq[query_loc]];
        int don = 0;

        //printf("%d %d Before_right\n", ref_loc, query_loc); 
        for(int i = 1; i < 100; i++){
            match_score = d_sub_mat[d_ref_seq[ref_loc+i]*5+d_query_seq[query_loc+i]];
            ext_score = ext_score + match_score;
            max_score = (ext_score > max_score) ? ext_score : max_score;

            if((max_score - ext_score) > y_drop && don == 0){
                don = 1;
                right_ext_score = ext_score;
            }
        }

        //printf("%d %d Before_left\n", ref_loc, query_loc); 
        ext_score = 0;
        max_score = 0;
        ref_loc = r_starts[id+thread_id] - 1;
        query_loc = q_starts[id+thread_id] - 1;
        ext_score = d_sub_mat[d_ref_seq[ref_loc]*5+d_query_seq[query_loc]];
        don = 0;
        
        for(int i = 1; i < 100; i++){
            match_score = d_sub_mat[d_ref_seq[ref_loc-i]*5+d_query_seq[query_loc-i]];
            ext_score = ext_score + match_score;
            max_score = (ext_score > max_score) ? ext_score : max_score;

            if((max_score - ext_score) > y_drop && don == 0){
                don = 1;
                left_ext_score = ext_score;
            }
        }

        total_score = seed_score + right_ext_score + left_ext_score; 
        if(total_score >= 3000){
            //printf("%d\n", total_score);
        }
        __syncthreads();
    }
}

int SeedAndFilter (std::vector<uint64_t> seed_offset_vector, int query_start){
    int ret = 0;
    cudaError_t err;

    uint64_t h_seed_offsets[13*cfg.chunk_size];
    
    uint32_t num_seeds = seed_offset_vector.size();
    assert(num_seeds <= 13*cfg.chunk_size);
    if (num_seeds == 0) {
        return ret;
    }

    gpu_lock.lock();
    for (uint32_t i = 0; i < num_seeds; i++) {
        h_seed_offsets[i] = seed_offset_vector[i];
    }
    
    err = cudaMemcpy(d_seed_offsets, h_seed_offsets, num_seeds*sizeof(uint64_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed!\n");
        exit(1);
    }

    int blockSize = 32;
    int numBlocks = 10;

    printf("Start find_hit_offsets\n");
    find_hit_offsets <<<numBlocks, blockSize>>> (num_seeds, query_start, d_ref_seq, d_query_seq, d_index_table, d_pos_table, d_seed_offsets, d_sub_mat);
    gpu_lock.unlock();

    return ret;
}

size_t InitializeProcessor (int t, int f){
    size_t ret = 0;
    cudaError_t err;

    err = cudaMalloc(&d_seed_offsets, 13*cfg.chunk_size*sizeof(uint64_t)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMalloc failed!\n");
        exit(1);
    }

    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMalloc failed!\n");
        exit(1);
    }
    
    int *sub_mat = (int *)malloc(25 * sizeof(int)); 

    sub_mat[24] = cfg.gact_sub_mat[10];
    for(int i = 0; i < 4; i++){
        sub_mat[i*5+4] = cfg.gact_sub_mat[10];
        sub_mat[4*5+i] = cfg.gact_sub_mat[10];
        sub_mat[i*5+i] = cfg.gact_sub_mat[i*4 + i - mat_offset[i]];
    }

    for(int i = 0; i < 4; i++){
        for(int j = i+1; j < 4; j++){
            sub_mat[i*5+j] = cfg.gact_sub_mat[i*4 + j - mat_offset[i]];
            sub_mat[j*5+i] = sub_mat[i*5+j];
        }
    }

    err = cudaMalloc(&d_sub_mat, 25*sizeof(int)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMalloc failed!\n");
        exit(1);
    }

    err = cudaMemcpy(d_sub_mat, sub_mat, 25*sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "1. Error: cudaMemcpy failed!\n");
        exit(1);
    }

    return ret;
}

void SendRefWriteRequest (size_t start_addr, size_t len){
    cudaError_t err;
    
    char* d_ref_seq_tmp;
    err = cudaMalloc(&d_ref_seq_tmp, len*sizeof(char)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMalloc failed!\n");
        exit(1);
    }

    err = cudaMemcpy(d_ref_seq_tmp, g_DRAM->buffer + start_addr, len*sizeof(char), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "1. Error: cudaMemcpy failed!\n");
        exit(1);
    }

    err = cudaMalloc(&d_ref_seq, len*sizeof(char)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMalloc failed!\n");
        exit(1);
    }
    
    int blockSize = 32;
    int numBlocks = (1 << 15);
    
    compress_string <<<numBlocks, blockSize>>> (len, d_ref_seq_tmp, d_ref_seq);
    
    uint32_t h_r[32];
    err = cudaMemcpy(h_r, d_ref_seq, 32, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "3. Error: cudaMemcpy failed!\n");
        exit(1);
    }
    
    for (int i = 0; i < 4; i++) {
        fprintf(stderr, "%d: 0x%x, ", i, h_r[i]);
    }

    fprintf(stderr, "\n");

    cudaFree(d_ref_seq_tmp);
}

void SendQueryWriteRequest (size_t start_addr, size_t len){
    cudaError_t err;
    
    char* d_query_seq_tmp;
    err = cudaMalloc(&d_query_seq_tmp, len*sizeof(char)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMemcpy(d_query_seq_tmp, g_DRAM->buffer + start_addr, len*sizeof(char), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed!\n");
        exit(1);
    }
    
    err = cudaMalloc(&d_query_seq, len*sizeof(char)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMalloc failed!\n");
        exit(1);
    }

    int blockSize = 32;
    int numBlocks = (1 << 15);
    
    compress_string <<<numBlocks, blockSize>>> (len, d_query_seq_tmp, d_query_seq);
    
    cudaFree(d_query_seq_tmp);
}

void SendSeedPosTable (uint32_t* index_table, uint32_t index_table_size, uint64_t* pos_table, uint32_t num_index){
    cudaError_t err;

    err = cudaMalloc(&d_index_table, index_table_size*sizeof(uint32_t)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMemcpy(d_index_table, index_table, index_table_size*sizeof(uint32_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_pos_table, num_index*sizeof(uint64_t)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMemcpy(d_pos_table, pos_table, num_index*sizeof(uint64_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed!\n");
        exit(1);
    }
    fprintf(stdout, "Sending seed position table to GPU successful.\n");
}

std::vector<tile_output> SendBatchRequest (std::vector<filter_tile> tiles, uint8_t align_fields, int thresh) {
        std::vector<tile_output> filtered_op;

            return filtered_op;
}

extend_output GACTXRequest (extend_tile tile, uint8_t align_fields) {
        extend_output op;

            return op;
}

void SendRequest (size_t ref_offset, size_t query_offset, size_t ref_length, size_t query_length, uint8_t align_fields){

}

void ShutdownProcessor(){
    cudaFree(d_ref_seq);
    cudaFree(d_query_seq);
    cudaFree(d_seed_offsets);
    cudaFree(d_index_table);
    cudaFree(d_pos_table);
}

DRAM *g_DRAM = nullptr;

InitializeProcessor_ptr g_InitializeProcessor = InitializeProcessor;
ShutdownProcessor_ptr g_ShutdownProcessor = ShutdownProcessor;
SendRequest_ptr g_SendRequest = SendRequest;
SendSeedPosTable_ptr g_SendSeedPosTable = SendSeedPosTable;
SeedAndFilter_ptr g_SeedAndFilter = SeedAndFilter;
SendRefWriteRequest_ptr g_SendRefWriteRequest = SendRefWriteRequest;
SendQueryWriteRequest_ptr g_SendQueryWriteRequest = SendQueryWriteRequest;       
SendBatchRequest_ptr g_SendBatchRequest = SendBatchRequest;
GACTXRequest_ptr g_GACTXRequest = GACTXRequest;
