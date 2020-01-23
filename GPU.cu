#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include "graph.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <algorithm>
#include <mutex>
#include <cstring>
#include <cstdlib>
#include "tbb/scalable_allocator.h"
#include <stdio.h>
#include <stdlib.h>

#define MAX_BLOCKS 1<<10
#define MAX_THREADS 1024 
#define TILE_SIZE 32
#define NUM_TILES 1000
#define MAX_HITS 7000000 
#define BLOCK_SIZE 128 
#define NUM_WARPS 4

std::mutex gpu_lock;

const char A_NT = 0;
const char C_NT = 1;
const char G_NT = 2;
const char T_NT = 3;
const char N_NT = 4;

int mat_offset[] = {0, 1, 3, 6};                                                             
int *d_sub_mat;

int err;                            
int check_status = 0;

int ref_len;
int query_len;
int seed_size;

struct timeval time1, time2, time3, time4, time5; 
long useconds1, seconds1, mseconds1;

char* d_ref_seq;
char* d_query_seq;
char* d_query_seq1;
char* d_query_rc_seq;
uint64_t* d_seed_offsets;
uint32_t* d_index_table;
uint32_t* d_num_seed_hits;
uint64_t* d_pos_table;
uint64_t* h_seed_offsets;
uint32_t* d_r_starts;
uint32_t* d_q_starts;
uint32_t* d_len;
long int* d_diag;
bool* d_done;
int *sub_mat;

uint32_t* h_r_starts;
uint32_t* h_q_starts;
int32_t* h_len;
bool* h_done;

__global__
void rev_comp (uint32_t len, char* src_seq, char* dst_seq){ 
    int thread_id = threadIdx.x;
    int block_dim = blockDim.x;
    int grid_dim = gridDim.x;
    int block_id = blockIdx.x;

    int stride = block_dim * grid_dim;
    uint32_t start = block_dim * block_id + thread_id;

    for (uint32_t i = start; i < len; i += stride) {
        char ch = src_seq[len - 1 -i];
        char dst = N_NT;
        if ((ch == 'a') || (ch == 'A'))
            dst = T_NT;
        else if ((ch == 'c') || (ch == 'C'))
            dst = G_NT;
        else if ((ch == 'g') || (ch == 'G'))
            dst = C_NT;
        else if ((ch == 't') || (ch == 'T'))
            dst = A_NT;
        dst_seq[i] = dst;
    }
}

__global__
void compress_string (uint32_t len, char* src_seq, char* dst_seq){ 
    int thread_id = threadIdx.x;
    int block_dim = blockDim.x;
    int grid_dim = gridDim.x;
    int block_id = blockIdx.x;

    int stride = block_dim * grid_dim;
    uint32_t start = block_dim * block_id + thread_id;

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
void find_num_hits (int num_seeds, const uint32_t* __restrict__ d_index_table, uint64_t* seed_offsets, uint32_t* d_num_seed_hits, int* seed_hit_num){

    int thread_id = threadIdx.x;
    int block_dim = blockDim.x;
    int grid_dim = gridDim.x;
    int block_id = blockIdx.x;

    int stride = block_dim * grid_dim;
    uint32_t start = block_dim * block_id + thread_id;

    uint32_t num_seed_hit;
    uint32_t seed;
    
    for (uint32_t id = start; id < num_seeds; id += stride) {
        seed = (seed_offsets[id] >> 32);

        // start and end from the seed block_id table
        num_seed_hit = d_index_table[seed];
        if (seed > 0){
            num_seed_hit -= d_index_table[seed-1];
        }

        d_num_seed_hits[id] = num_seed_hit;
        seed_hit_num[id] = num_seed_hit;
    }
}

__global__
void find_anchors2 (int num_seeds, const char* __restrict__  d_ref_seq, const char* __restrict__  d_query_seq, const uint32_t* __restrict__  d_index_table, const uint64_t* __restrict__ d_pos_table, uint64_t*  d_seed_offsets, int *d_sub_mat, int xdrop, int xdrop_threshold, uint32_t* d_r_starts, uint32_t* d_q_starts, uint32_t* d_len, long int* d_diag, bool* d_done, int ref_len, int query_len, int seed_size, int* seed_hit_num, int num_hits, bool rev){

    int thread_id = threadIdx.x;
    int block_id = blockIdx.x;
    int warp_size = warpSize;
    int lane_id = thread_id%warp_size;
    int warp_id = (thread_id-lane_id)/warp_size;

    __shared__ uint32_t start, end;
    __shared__ uint32_t seed;
    __shared__ uint64_t seed_offset;

    __shared__ uint32_t ref_loc[NUM_WARPS];
    __shared__ uint32_t query_loc;
    __shared__ int total_score[NUM_WARPS];
    __shared__ int prev_score[NUM_WARPS];
    __shared__ int prev_max_score[NUM_WARPS];
    __shared__ bool right_edge[NUM_WARPS]; 
    __shared__ bool left_edge[NUM_WARPS]; 
    __shared__ bool right_xdrop_found[NUM_WARPS]; 
    __shared__ bool left_xdrop_found[NUM_WARPS]; 
    __shared__ uint32_t left_extent[NUM_WARPS];
    __shared__ uint32_t right_extent[NUM_WARPS];
    __shared__ uint32_t tile[NUM_WARPS];
    __shared__ char seed_query[32];

    int thread_score;
    int max_thread_score;
    bool xdrop_done;
    int temp;
    int ref_pos;
    int query_pos;

    __shared__ int sub_mat[25];

    if(thread_id < 25){
        sub_mat[thread_id] = d_sub_mat[thread_id];
    }

    if(thread_id == 0){
        seed_offset = d_seed_offsets[block_id];
        seed = (seed_offset >> 32);
        query_loc = ((seed_offset << 32) >> 32);

        // start and end from the seed block_id table
        end = d_index_table[seed];
        start = 0;
        if (seed > 0){
            start = d_index_table[seed-1];
        }
    }

    if(warp_id == 0)
        seed_query[lane_id] = d_query_seq[query_loc+lane_id]; 
    __syncthreads();

    for (int id1 = start; id1 < end; id1 += NUM_WARPS) {
        if(id1+warp_id < end){ 
            if(lane_id == 0){ 
                ref_loc[warp_id]   = d_pos_table[id1+warp_id];
                total_score[warp_id] = 0; 
                tile[warp_id] = 0;
            }

            //////////////////////////////////////////////////////////////////

            thread_score = 0;
            if(lane_id < seed_size){
                thread_score = sub_mat[d_ref_seq[ref_loc[warp_id]+lane_id]*5+seed_query[lane_id]];
            }

            for (int offset = warp_size >> 1; offset > 0; offset = offset >> 1)
                thread_score += __shfl_down_sync(0x13, thread_score, offset);

            if(lane_id == 0){
                total_score[warp_id] += thread_score;
            }
            __syncwarp();

            ////////////////////////////////////////////////////////////////

            tile[warp_id] = 0;
            right_xdrop_found[warp_id] = false;
            right_edge[warp_id] = false;
            prev_score[warp_id] = 0;
            prev_max_score[warp_id] = 0;
            right_extent[warp_id] = 0;

            while(!right_xdrop_found[warp_id] && !right_edge[warp_id]){
                ref_pos   = ref_loc[warp_id]   + seed_size + lane_id + tile[warp_id]*warp_size;
                query_pos = query_loc + seed_size + lane_id + tile[warp_id]*warp_size;

                if(ref_pos < ref_len && query_pos < query_len){
                    thread_score = sub_mat[d_ref_seq[ref_pos]*5+d_query_seq[query_pos]];
                }

#pragma unroll
                for (int offset = 1; offset < warp_size; offset = offset << 1){
                    temp = __shfl_up_sync(0xFFFFFFFF, thread_score, offset);

                    if(lane_id >= offset){
                        thread_score += temp;
                    }
                }

                thread_score += prev_score[warp_id];
                max_thread_score = max(thread_score, prev_max_score[warp_id]);
                __syncwarp();

#pragma unroll
                for (int offset = 1; offset < warp_size; offset = offset << 1){
                    temp = __shfl_up_sync(0xFFFFFFFF, max_thread_score, offset);

                    if(lane_id >= offset){
                        max_thread_score = max(max_thread_score, temp);
                    }
                }

                xdrop_done = ((max_thread_score-thread_score) > xdrop);
                __syncwarp();

#pragma unroll
                for (int offset = 1; offset < warp_size; offset = offset << 1){
                    xdrop_done |= __shfl_up_sync(0xFFFFFFFF, xdrop_done, offset);
                }

                if(lane_id == warp_size-1){
                    if(xdrop_done){
                        total_score[warp_id]+=max_thread_score;
                        right_xdrop_found[warp_id] = true;
                        right_extent[warp_id] = tile[warp_id]*warp_size-1;
                    }
                    else if(ref_pos > ref_len || query_pos > query_len)
                        right_edge[warp_id] = true;
                    else{
                        prev_score[warp_id] = thread_score;
                        prev_max_score[warp_id] = max_thread_score;
                    }
                }
                __syncwarp();

                tile[warp_id]++;
            }

            ////////////////////////////////////////////////////////////////

            tile[warp_id] = 0;
            left_xdrop_found[warp_id] = false;
            left_edge[warp_id] = false;
            prev_score[warp_id] = 0;
            prev_max_score[warp_id] = 0;
            left_extent[warp_id] = 0;

            while(!left_xdrop_found[warp_id] && !left_edge[warp_id]){

                ref_pos   = ref_loc[warp_id] - lane_id - 1 - tile[warp_id]*warp_size;
                query_pos = query_loc - lane_id - 1 - tile[warp_id]*warp_size;

                if(ref_pos >= 0  && query_pos >= 0){
                    thread_score = sub_mat[d_ref_seq[ref_pos]*5+d_query_seq[query_pos]];
                }

#pragma unroll
                for (int offset = 1; offset < warp_size; offset = offset << 1){
                    temp = __shfl_up_sync(0xFFFFFFFF, thread_score, offset);

                    if(lane_id >= offset){
                        thread_score += temp;
                    }
                }

                thread_score += prev_score[warp_id];
                max_thread_score = max(thread_score, prev_max_score[warp_id]);
                __syncwarp();

#pragma unroll
                for (int offset = 1; offset < warp_size; offset = offset << 1){
                    temp = __shfl_up_sync(0xFFFFFFFF, max_thread_score, offset);

                    if(lane_id >= offset){
                        max_thread_score = max(max_thread_score, temp);
                    }
                }

                xdrop_done = ((max_thread_score-thread_score) > xdrop);
                __syncwarp();

#pragma unroll
                for (int offset = 1; offset < warp_size; offset = offset << 1){
                    xdrop_done |= __shfl_up_sync(0xFFFFFFFF, xdrop_done, offset);
                }

                if(lane_id == warp_size-1){
                    if(xdrop_done){
                        total_score[warp_id]+=max_thread_score;
                        left_xdrop_found[warp_id] = true;
                        left_extent[warp_id] = tile[warp_id]*warp_size-1;
                    }
                    else if(ref_pos < 0 || query_pos < 0)
                        left_edge[warp_id] = true;
                    else{
                        prev_score[warp_id] = thread_score;
                        prev_max_score[warp_id] = max_thread_score;
                    }
                }
                __syncwarp();

                tile[warp_id]++;
            }

            //////////////////////////////////////////////////////////////////

            if(lane_id == 0){
                int dram_address = seed_hit_num[block_id]-id1 - warp_id+start-1;

                if(total_score[warp_id] >= xdrop_threshold){
                    d_r_starts[dram_address] = ref_loc[warp_id] - left_extent[warp_id];
                    d_q_starts[dram_address] = query_loc - left_extent[warp_id];
                    d_len[dram_address] = left_extent[warp_id]+right_extent[warp_id]+seed_size;
                    d_diag[dram_address] = ref_loc[warp_id] - query_loc;
                    d_done[dram_address] = true;
                }
                else{
                    d_r_starts[dram_address] = 0;
                    d_q_starts[dram_address] = 0;
                    d_len[dram_address] = 0;
                    d_diag[dram_address] = 0;
                    d_done[dram_address] = false;
                }
            }
            __syncwarp();
        }
    }
}

__global__
void find_anchors0 (int num_seeds, const char* __restrict__  d_ref_seq, const char* __restrict__  d_query_seq, const uint32_t* __restrict__  d_index_table, const uint64_t* __restrict__ d_pos_table, uint64_t*  d_seed_offsets, int *d_sub_mat, int xdrop, int xdrop_threshold, uint32_t* d_r_starts, uint32_t* d_q_starts, uint32_t* d_len, long int* d_diag,bool* d_done, int ref_len, int query_len, int seed_size, int* seed_hit_num, int num_hits, bool rev){

    int thread_id = threadIdx.x;
    int block_id = blockIdx.x;
    int warp_size = warpSize;
    int lane_id = thread_id%warp_size;
    int warp_id = (thread_id-lane_id)/warp_size;

    __shared__ uint32_t start, end;
    __shared__ uint32_t seed;
    __shared__ uint64_t seed_offset;

    __shared__ uint32_t ref_loc[NUM_WARPS];
    __shared__ uint32_t query_loc;
    __shared__ int total_score[NUM_WARPS];
    __shared__ int prev_score[NUM_WARPS];
    __shared__ int prev_max_score[NUM_WARPS];
    __shared__ bool right_edge[NUM_WARPS]; 
    __shared__ bool left_edge[NUM_WARPS]; 
    __shared__ bool right_xdrop_found[NUM_WARPS]; 
    __shared__ bool left_xdrop_found[NUM_WARPS]; 
    __shared__ uint32_t left_extent[NUM_WARPS];
    __shared__ uint32_t right_extent[NUM_WARPS];
    __shared__ uint32_t tile[NUM_WARPS];

    int thread_score;
    int max_thread_score;
    bool xdrop_done;
    int temp;
    int ref_pos;
    int query_pos;

    __shared__ int sub_mat[25];

    if(thread_id < 25){
        sub_mat[thread_id] = d_sub_mat[thread_id];
    }

    if(thread_id == 0){
        seed_offset = d_seed_offsets[block_id];
        seed = (seed_offset >> 32);
        query_loc = ((seed_offset << 32) >> 32);

        // start and end from the seed block_id table
        end = d_index_table[seed];
        start = 0;
        if (seed > 0){
            start = d_index_table[seed-1];
        }
    }
    __syncthreads();

    for (int id1 = start; id1 < end; id1 += NUM_WARPS) {
        if(id1+warp_id < end){ 
            if(lane_id == 0){ 
                ref_loc[warp_id]   = d_pos_table[id1+warp_id];
                total_score[warp_id] = 0; 
                tile[warp_id] = 0;
            }

            //////////////////////////////////////////////////////////////////

            thread_score = 0;
            if(lane_id < seed_size){
                thread_score = sub_mat[d_ref_seq[ref_loc[warp_id]+lane_id]*5+d_query_seq[query_loc+lane_id]];
            }

            for (int offset = warp_size >> 1; offset > 0; offset = offset >> 1)
                thread_score += __shfl_down_sync(0x13, thread_score, offset);

            if(lane_id == 0){
                total_score[warp_id] += thread_score;
            }
            __syncwarp();

            ////////////////////////////////////////////////////////////////

            tile[warp_id] = 0;
            right_xdrop_found[warp_id] = false;
            right_edge[warp_id] = false;
            prev_score[warp_id] = 0;
            prev_max_score[warp_id] = 0;
            right_extent[warp_id] = 0;

            while(tile[warp_id] < NUM_TILES && !right_xdrop_found[warp_id] && !right_edge[warp_id]){
                ref_pos   = ref_loc[warp_id]   + seed_size + lane_id + tile[warp_id]*warp_size;
                query_pos = query_loc + seed_size + lane_id + tile[warp_id]*warp_size;

                if(ref_pos < ref_len && query_pos < query_len){
                    thread_score = sub_mat[d_ref_seq[ref_pos]*5+d_query_seq[query_pos]];
                }

#pragma unroll
                for (int offset = 1; offset < warp_size; offset = offset << 1){
                    temp = __shfl_up_sync(0xFFFFFFFF, thread_score, offset);

                    if(lane_id >= offset){
                        thread_score += temp;
                    }
                }

                thread_score += prev_score[warp_id];
                max_thread_score = max(thread_score, prev_max_score[warp_id]);
                __syncwarp();

#pragma unroll
                for (int offset = 1; offset < warp_size; offset = offset << 1){
                    temp = __shfl_up_sync(0xFFFFFFFF, max_thread_score, offset);

                    if(lane_id >= offset){
                        max_thread_score = max(max_thread_score, temp);
                    }
                }

                xdrop_done = ((max_thread_score-thread_score) > xdrop);
                __syncwarp();

#pragma unroll
                for (int offset = 1; offset < warp_size; offset = offset << 1){
                    xdrop_done |= __shfl_up_sync(0xFFFFFFFF, xdrop_done, offset);
                }

                if(lane_id == warp_size-1){
                    if(xdrop_done){
                        total_score[warp_id]+=max_thread_score;
                        right_xdrop_found[warp_id] = true;
                        right_extent[warp_id] = tile[warp_id]*TILE_SIZE-1;
                    }
                    else if(ref_pos > ref_len || query_pos > query_len)
                        right_edge[warp_id] = true;
                    else{
                        prev_score[warp_id] = thread_score;
                        prev_max_score[warp_id] = max_thread_score;
                    }
                }
                __syncwarp();

                tile[warp_id]++;
            }

            ////////////////////////////////////////////////////////////////

            tile[warp_id] = 0;
            left_xdrop_found[warp_id] = false;
            left_edge[warp_id] = false;
            left_extent[warp_id] = 0;
            prev_score[warp_id] = 0;
            prev_max_score[warp_id] = 0;

            while(tile[warp_id] < NUM_TILES && !left_xdrop_found[warp_id] && !left_edge[warp_id]){

                ref_pos   = ref_loc[warp_id] - lane_id - 1 - tile[warp_id]*warp_size;
                query_pos = query_loc - lane_id - 1 - tile[warp_id]*warp_size;

                if(ref_pos >= 0  && query_pos >= 0){
                    thread_score = sub_mat[d_ref_seq[ref_pos]*5+d_query_seq[query_pos]];
                }

#pragma unroll
                for (int offset = 1; offset < warp_size; offset = offset << 1){
                    temp = __shfl_up_sync(0xFFFFFFFF, thread_score, offset);

                    if(lane_id >= offset){
                        thread_score += temp;
                    }
                }

                thread_score += prev_score[warp_id];
                max_thread_score = max(thread_score, prev_max_score[warp_id]);
                __syncwarp();

#pragma unroll
                for (int offset = 1; offset < warp_size; offset = offset << 1){
                    temp = __shfl_up_sync(0xFFFFFFFF, max_thread_score, offset);

                    if(lane_id >= offset){
                        max_thread_score = max(max_thread_score, temp);
                    }
                }

                xdrop_done = ((max_thread_score-thread_score) > xdrop);
                __syncwarp();

#pragma unroll
                for (int offset = 1; offset < warp_size; offset = offset << 1){
                    xdrop_done |= __shfl_up_sync(0xFFFFFFFF, xdrop_done, offset);
                }

                if(lane_id == warp_size-1){
                    if(xdrop_done){
                        total_score[warp_id]+=max_thread_score;
                        left_xdrop_found[warp_id] = true;
                        left_extent[warp_id] = tile[warp_id]*BLOCK_SIZE-1;
                    }
                    else if(ref_pos < 0 || query_pos < 0)
                        left_edge[warp_id] = true;
                    else{
                        prev_score[warp_id] = thread_score;
                        prev_max_score[warp_id] = max_thread_score;
                    }
                }
                __syncwarp();

                tile[warp_id]++;
            }

            //////////////////////////////////////////////////////////////////

            if(lane_id == 0){
                int dram_address = seed_hit_num[block_id]-id1 - warp_id+start-1;

                if(total_score[warp_id] < 0)
                    printf("t %d %d\n", total_score[warp_id], rev);

                if(!left_edge[warp_id] && !right_edge[warp_id]){
                    if(right_xdrop_found[warp_id] && left_xdrop_found[warp_id]){
                        if(total_score[warp_id] >= xdrop_threshold){
                            d_r_starts[dram_address] = ref_loc[warp_id] - left_extent[warp_id];
                            d_q_starts[dram_address] = query_loc - left_extent[warp_id];
                            d_len[dram_address] = left_extent[warp_id]+right_extent[warp_id]+seed_size;
                            d_done[dram_address] = true;
                        }
                        else{
                            d_r_starts[dram_address] = 0;
                            d_q_starts[dram_address] = 0;
                            d_len[dram_address] = 1;
                            d_done[dram_address] = false;
                        }
                    }
                    else{
                        d_r_starts[dram_address] = ref_loc[warp_id];
                        d_q_starts[dram_address] = query_loc;
                        d_len[dram_address] = seed_size;
                        d_done[dram_address] = false;
                    }
                }
                else{
                    d_r_starts[dram_address] = 0;
                    d_q_starts[dram_address] = 0;
                    d_len[dram_address] = 2;
                    d_done[dram_address] = false;
                }
            }
            __syncwarp();
        }
    }
}

int SeedAndFilter (std::vector<uint64_t> seed_offset_vector, bool rev){

    gpu_lock.lock();
    int ret = 0;
    cudaError_t err;
    seed_size = 19;

    uint32_t num_hits;
    int total_anchors = 0;
    int total1 = 0;

    uint32_t num_seeds = seed_offset_vector.size();
    assert(num_seeds <= 13*cfg.chunk_size);

    if (num_seeds == 0) {
        return ret;
    }

    for (uint32_t i = 0; i < num_seeds; i++) {
        h_seed_offsets[i] = seed_offset_vector[i];
    }

    err = cudaMemcpy(d_seed_offsets, h_seed_offsets, num_seeds*sizeof(uint64_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed!\n");
        exit(1);
    }

    thrust::device_vector<int> seed_hit_num(num_seeds);
    int* seed_hit_num_array = thrust::raw_pointer_cast(&seed_hit_num[0]);

    find_num_hits <<<MAX_BLOCKS, MAX_THREADS>>> (num_seeds, d_index_table, d_seed_offsets, d_num_seed_hits, seed_hit_num_array);

    thrust::inclusive_scan(seed_hit_num.begin(), seed_hit_num.end(), seed_hit_num.begin());

    err = cudaMemcpy(&num_hits, (seed_hit_num_array+num_seeds-1), sizeof(uint32_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed! SF1\n");
        exit(1);
    }

    gettimeofday(&time2, NULL);

    if(rev)
        find_anchors2 <<<num_seeds,BLOCK_SIZE>>> (num_seeds, d_ref_seq, d_query_rc_seq, d_index_table, d_pos_table, d_seed_offsets, d_sub_mat, cfg.xdrop, cfg.xdrop_threshold, d_r_starts, d_q_starts, d_len, d_diag, d_done, ref_len, query_len, seed_size, seed_hit_num_array, num_hits, rev);
    else
        find_anchors2 <<<num_seeds,BLOCK_SIZE>>> (num_seeds, d_ref_seq, d_query_seq, d_index_table, d_pos_table, d_seed_offsets, d_sub_mat, cfg.xdrop, cfg.xdrop_threshold, d_r_starts, d_q_starts, d_len, d_diag, d_done, ref_len, query_len, seed_size, seed_hit_num_array, num_hits, rev);

//    if(rev)
//        find_anchors0 <<<num_seeds,BLOCK_SIZE>>> (num_seeds, d_ref_seq, d_query_rc_seq, d_index_table, d_pos_table, d_seed_offsets, d_sub_mat, cfg.xdrop, cfg.xdrop_threshold, d_r_starts, d_q_starts, d_len, d_done, ref_len, query_len, seed_size, seed_hit_num_array, num_hits, rev);
//    else
//        find_anchors0 <<<num_seeds,BLOCK_SIZE>>> (num_seeds, d_ref_seq, d_query_seq, d_index_table, d_pos_table, d_seed_offsets, d_sub_mat, cfg.xdrop, cfg.xdrop_threshold, d_r_starts, d_q_starts, d_len, d_done, ref_len, query_len, seed_size, seed_hit_num_array, num_hits, rev);

    err = cudaMemcpy(h_r_starts, d_r_starts, num_hits*sizeof(uint32_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed! SF6\n");
        exit(1);
    }
    
    err = cudaMemcpy(h_q_starts, d_q_starts, num_hits*sizeof(uint32_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed! SF7\n");
        exit(1);
    }
    
    err = cudaMemcpy(h_len, d_len, num_hits*sizeof(uint32_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed! SF8\n");
        exit(1);
    }

    err = cudaMemcpy(h_done, d_done, num_hits*sizeof(bool), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed! SF9\n");
        exit(1);
    }

    gpu_lock.unlock();

    gettimeofday(&time4, NULL);

    for(int i = 0; i < num_hits; i++){
        if(h_done[i]){
            fprintf(stdout, "%d %d\n", h_r_starts[i], h_q_starts[i]);
        }

        if(h_done[i]){
            total_anchors++;
        }
        else{
            total1++;
        }
    }

    useconds1 = time4.tv_usec - time2.tv_usec;
    seconds1  = time4.tv_sec  - time2.tv_sec;
    mseconds1 = ((seconds1) * 1000 + useconds1/1000.0) + 0.5;
//    fprintf(stdout, "%d %d %d %d\n", rev, total1, total_anchors, num_hits);

    return ret;
}

size_t InitializeProcessor (int t, int f){
    size_t ret = 0;
    cudaError_t err;

    h_r_starts = (uint32_t*) calloc(MAX_HITS, sizeof(uint32_t));
    h_q_starts = (uint32_t*) calloc(MAX_HITS, sizeof(uint32_t));
    h_len      = (int32_t*) calloc(MAX_HITS, sizeof(uint32_t));
    h_done     = (bool*) calloc(MAX_HITS, sizeof(bool));

    err = cudaMalloc(&d_r_starts, MAX_HITS*sizeof(uint32_t)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "2 Error: cudaMalloc failed! SF2\n");
        exit(1);
    }

    err = cudaMalloc(&d_q_starts, MAX_HITS*sizeof(uint32_t)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "3 Error: cudaMalloc failed! SF3\n");
        exit(1);
    }

    err = cudaMalloc(&d_len, MAX_HITS*sizeof(uint32_t)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "4 Error: cudaMalloc failed! SF4\n");
        exit(1);
    }

    err = cudaMalloc(&d_diag, MAX_HITS*sizeof(long int)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "4 Error: cudaMalloc failed! SF4\n");
        exit(1);
    }

    err = cudaMalloc(&d_done, MAX_HITS*sizeof(bool)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "5 Error: cudaMalloc failed! SF5\n");
        exit(1);
    }
    err = cudaMalloc(&d_seed_offsets, 13*cfg.chunk_size*sizeof(uint64_t)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "1 Error: cudaMalloc failed!\n");
        exit(1);
    }

    err = cudaMalloc(&d_num_seed_hits, 13*cfg.chunk_size*sizeof(uint32_t));
    if (err != cudaSuccess) {
        fprintf(stderr, "2 Error: cudaMalloc failed!\n");
        exit(1);
    }

    sub_mat = (int *)malloc(25 * sizeof(int)); 

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
    ref_len = len;
    
    char* d_ref_seq_tmp;
    err = cudaMalloc(&d_ref_seq_tmp, len*sizeof(char)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "1 Error: cudaMalloc failed!\n");
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
    
    compress_string <<<MAX_BLOCKS, MAX_THREADS>>> (len, d_ref_seq_tmp, d_ref_seq);
    
    cudaFree(d_ref_seq_tmp);
}

void SendQueryWriteRequest (size_t start_addr, size_t len){
    cudaError_t err;
    query_len = len;
    
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

    err = cudaMalloc(&d_query_rc_seq, len*sizeof(char)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMalloc failed!\n");
        exit(1);
    }

    rev_comp <<<MAX_BLOCKS, MAX_THREADS>>> (len, d_query_seq_tmp, d_query_rc_seq);
    compress_string <<<MAX_BLOCKS, MAX_THREADS>>> (len, d_query_seq_tmp, d_query_seq);

    cudaFree(d_query_seq_tmp);
}

void SendSeedPosTable (uint32_t* index_table, uint32_t index_table_size, uint64_t* pos_table, uint32_t num_index){
    cudaError_t err;

    h_seed_offsets = (uint64_t*) malloc(13*cfg.chunk_size*sizeof(uint64_t));

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
    cudaFree(d_query_rc_seq);
    cudaFree(d_seed_offsets);
    cudaFree(d_index_table);
    cudaFree(d_pos_table);
    cudaFree(d_r_starts);
    cudaFree(d_q_starts);
    cudaFree(d_len);
    cudaFree(d_done);
    free(h_r_starts);
    free(h_q_starts);
    free(h_len);
    free(h_done);
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
