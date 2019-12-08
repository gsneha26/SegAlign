#include "graph.h"
#include <algorithm>
#include <mutex>
#include <cstring>
#include <cstdlib>
#include "tbb/scalable_allocator.h"
#include <stdio.h>
#include <stdlib.h>

#define NUM_THREADS 64 
#define NUM_BLOCKS 2048 
#define HIT_LIMIT 512
#define OUT_LIMIT 1024

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

char* d_ref_seq;
char* d_query_seq;
uint64_t* d_seed_offsets;
uint32_t* d_index_table;
uint64_t* d_pos_table;
uint32_t* d_r_starts;
uint32_t* d_q_starts;
uint32_t* d_len;
bool* d_done;
uint64_t* h_seed_offsets;
int *sub_mat;

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
void find_anchors2 (int num_seeds, char* d_ref_seq, char* d_query_seq, uint32_t* d_index_table, uint64_t* d_pos_table, uint64_t* d_seed_offsets, int *d_sub_mat, int xdrop, int xdrop_threshold, uint32_t* d_r_starts, uint32_t* d_q_starts, uint32_t* d_len, bool* d_done, int ref_len, int query_len, int seed_size){

    int thread_id = threadIdx.x;
    int block_id = blockIdx.x;

    __shared__ uint32_t start, end;
    __shared__ uint32_t q_start;
    __shared__ uint32_t seed;
    __shared__ uint64_t seed_offset;

    __shared__ uint32_t ref_loc;
    __shared__ uint32_t query_loc;
    __shared__ int total_score;
    __shared__ int warp_size;
    int max_score;
    __shared__ bool right_xdrop_found; 
    __shared__ bool left_xdrop_found; 
    __shared__ uint32_t left_extent;
    __shared__ uint32_t right_extent;
    int thread_score;

    __shared__ int score[NUM_THREADS];
    __shared__ int sub_mat[25];
    __shared__ uint32_t total_anchors;

    if(thread_id < 25){
        sub_mat[thread_id] = d_sub_mat[thread_id];
    }

    if(thread_id == 0){
        total_anchors = 0;
        warp_size = warpSize;
        seed_offset = d_seed_offsets[block_id];
    }

    int lane_id = thread_id % warp_size;
    __syncthreads();

    seed = (seed_offset >> 32);
    q_start = ((seed_offset << 32) >> 32);

    // start and end from the seed block_id table
    end = d_index_table[seed];
    start = 0;
    if (seed > 0){
        start = d_index_table[seed-1];
    }
    else{
        start = 0;
        end = 0;
    }

    for (int id1 = start; id1 < end; id1 += 1) {
        ref_loc   = d_pos_table[0];
        query_loc = q_start;
        total_score = 0; 

        //////////////////////////////////////////////////////////////////

        thread_score = 0;
        if(thread_id < seed_size){
            thread_score = sub_mat[d_ref_seq[ref_loc+thread_id]*5+d_query_seq[query_loc+thread_id]];
        }

        for (int offset = warp_size/2; offset > 0; offset /= 2)
            thread_score += __shfl_down_sync(0x13, thread_score, offset);

        if(thread_id == 0){
            total_score += thread_score;
        }

        //////////////////////////////////////////////////////////////////

        thread_score = 0;
        if(ref_loc+seed_size+thread_id < ref_len && query_loc+seed_size+thread_id < query_len){
            score[thread_id] = sub_mat[d_ref_seq[ref_loc+seed_size+thread_id]*5+d_query_seq[query_loc+seed_size+thread_id]];
        }
        __syncthreads();

        //parallel prefix sum
        int k_val;
        for(int k = 1; k < NUM_THREADS; k=k*2){
            if(thread_id >= k){
                k_val = score[thread_id-k];
            }
            __syncthreads();
            if(thread_id >= k){
                score[thread_id] += k_val;
            }
            __syncthreads();
        }

        if(thread_id == 0){
            right_xdrop_found = false;
            max_score = 0;
            for(int i = 0; i < NUM_THREADS; i++){
                if(score[i] > max_score){
                    max_score = score[i];
                }

                if(max_score-score[i] > xdrop && right_xdrop_found == false){
                    total_score+=max_score;
                    right_xdrop_found = true;
                    right_extent = i;
                    break;
                }
            }
        }
        __syncthreads();

        //////////////////////////////////////////////////////////////////

        thread_score = 0;
        if(ref_loc >= thread_id+1 && query_loc>= thread_id+1){
            score[thread_id] = sub_mat[d_ref_seq[ref_loc-thread_id-1]*5+d_query_seq[query_loc-thread_id-1]];
        }
        __syncthreads();

        //parallel prefix sum
        for(int k = 1; k < NUM_THREADS; k=k*2){
            if(thread_id >= k){
                k_val = score[thread_id-k];
            }
            __syncthreads();
            if(thread_id >= k){
                score[thread_id] += k_val;
            }
            __syncthreads();
        }

        if(thread_id == 0){
            left_xdrop_found = false;
            max_score = 0;
            for(int i = 0; i < NUM_THREADS; i++){
                if(score[i] > max_score){
                    max_score = score[i];
                }

                if(max_score-score[i] > xdrop && left_xdrop_found == false){
                    total_score+=max_score;
                    left_xdrop_found = true;
                    left_extent = i;
                    break;
                }
            }
        }
        __syncthreads();

        //////////////////////////////////////////////////////////////////

        if(total_anchors < OUT_LIMIT){
            if(right_xdrop_found && left_xdrop_found){
                if(total_score >= xdrop_threshold){
//                    d_r_starts[block_id*OUT_LIMIT + total_anchors] = ref_loc - left_extent;
//                    d_q_starts[block_id*OUT_LIMIT + total_anchors] = query_loc - left_extent;
//                    d_len[block_id*OUT_LIMIT + total_anchors] = left_extent+right_extent+seed_size;
//                    d_done[block_id*OUT_LIMIT + total_anchors] = true;
                    d_r_starts[block_id] = ref_loc - left_extent;
                    d_q_starts[block_id] = query_loc - left_extent;
                    d_len[block_id] = left_extent+right_extent+seed_size;
                    d_done[block_id] = true;
                    total_anchors++;
                }
            }
            else{
                d_r_starts[block_id] = ref_loc;
                d_q_starts[block_id] = query_loc;
                d_len[block_id] = seed_size;
                d_done[block_id] = false;
                total_anchors++;
            }
        }

//        if(block_id == 0 && thread_id == 0){
//            printf("%d %d %d\n", start, total_score, score[0]);
//        }

        __syncthreads();
    }
}

__global__
void find_anchors1 (int num_seeds, char* d_ref_seq, char* d_query_seq, uint32_t* d_index_table, uint64_t* d_pos_table, uint64_t* d_seed_offsets, int *d_sub_mat, int xdrop, int xdrop_threshold, uint32_t* d_r_starts, uint32_t* d_q_starts, uint32_t* d_len, bool* d_done, int ref_len, int query_len, int seed_size){

    int thread_id = threadIdx.x;
    int block_id = blockIdx.x;

    __shared__ uint32_t start, end;
    __shared__ uint32_t q_start;
    __shared__ uint32_t seed;
    __shared__ uint64_t seed_offset;

    __shared__ uint32_t ref_loc;
    __shared__ uint32_t query_loc;
    __shared__ int total_score;
    __shared__ int warp_size;
    int max_score;
    __shared__ bool right_xdrop_found; 
    __shared__ bool left_xdrop_found; 
    __shared__ uint32_t left_extent;
    __shared__ uint32_t right_extent;
    int thread_score;
    int temp;

    __shared__ int score[NUM_THREADS];
    __shared__ int sub_mat[25];
    __shared__ uint32_t total_anchors;

    if(thread_id < 25){
        sub_mat[thread_id] = d_sub_mat[thread_id];
    }

    if(thread_id == 0){
        total_anchors = 0;
        warp_size = warpSize;
        seed_offset = d_seed_offsets[block_id];
    }

    int lane_id = thread_id % warp_size;
    __syncthreads();

    seed = (seed_offset >> 32);
    q_start = ((seed_offset << 32) >> 32);

    // start and end from the seed block_id table
    end = d_index_table[seed];
    start = 0;
    if (seed > 0){
        start = d_index_table[seed-1];
    }
    else{
        start = 0;
        end = 0;
    }

    for (int id1 = start; id1 < end; id1 += 1) {
        ref_loc   = d_pos_table[0];
        query_loc = q_start;
        total_score = 0; 

        //////////////////////////////////////////////////////////////////

        thread_score = 0;
        if(thread_id < seed_size){
            thread_score = sub_mat[d_ref_seq[ref_loc+thread_id]*5+d_query_seq[query_loc+thread_id]];
        }

        for (int offset = warp_size/2; offset > 0; offset /= 2)
            thread_score += __shfl_down_sync(0x13, thread_score, offset);

        if(thread_id == 0){
            total_score += thread_score;
        }

        //////////////////////////////////////////////////////////////////

        thread_score = 0;
        if(ref_loc+seed_size+thread_id < ref_len && query_loc+seed_size+thread_id < query_len){
            thread_score = sub_mat[d_ref_seq[ref_loc+seed_size+thread_id]*5+d_query_seq[query_loc+seed_size+thread_id]];
        }

        for (int offset = 1; offset < warp_size; offset *= 2){
            temp = __shfl_up_sync(0xFFFFFFFF, thread_score, offset);

            if(lane_id >= offset){
                thread_score += temp;
            }
        }
        __syncthreads();

        score[thread_id] = thread_score;
        __syncthreads();

        if(thread_id >= warp_size){
            score[thread_id] += score[warp_size-1];
        }
        __syncthreads();

        if(thread_id == 0){
            right_xdrop_found = false;
            max_score = 0;
            for(int i = 0; i < NUM_THREADS; i++){
                if(score[i] > max_score){
                    max_score = score[i];
                }

                if(max_score-score[i] > xdrop && right_xdrop_found == false){
                    total_score+=max_score;
                    right_xdrop_found = true;
                    right_extent = i;
                    break;
                }
            }
        }
        __syncthreads();

        //////////////////////////////////////////////////////////////////

        thread_score = 0;
        if(ref_loc >= thread_id+1 && query_loc>= thread_id+1){
            thread_score = sub_mat[d_ref_seq[ref_loc-thread_id-1]*5+d_query_seq[query_loc-thread_id-1]];
        }

        for (int offset = 1; offset < warp_size; offset *= 2){
            temp = __shfl_up_sync(0xFFFFFFFF, thread_score, offset);

            if(lane_id >= offset){
                thread_score += temp;
            }
        }

        score[thread_id] = thread_score;
        __syncthreads();

        if(thread_id >= warp_size){
            score[thread_id] += score[warp_size-1];
        }
        __syncthreads();

        if(thread_id == 0){
            left_xdrop_found = false;
            max_score = 0;
            for(int i = 0; i < NUM_THREADS; i++){
                if(score[i] > max_score){
                    max_score = score[i];
                }

                if(max_score-score[i] > xdrop && left_xdrop_found == false){
                    total_score+=max_score;
                    left_xdrop_found = true;
                    left_extent = i;
                    break;
                }
            }
        }
        __syncthreads();

        //////////////////////////////////////////////////////////////////

        if(total_anchors < OUT_LIMIT){
            if(right_xdrop_found && left_xdrop_found){
                if(total_score >= xdrop_threshold){
                    printf("here\n");
//                    d_r_starts[block_id*OUT_LIMIT + total_anchors] = ref_loc - left_extent;
//                    d_q_starts[block_id*OUT_LIMIT + total_anchors] = query_loc - left_extent;
//                    d_len[block_id*OUT_LIMIT + total_anchors] = left_extent+right_extent+seed_size;
//                    d_done[block_id*OUT_LIMIT + total_anchors] = true;
                    d_r_starts[block_id] = ref_loc - left_extent;
                    d_q_starts[block_id] = query_loc - left_extent;
                    d_len[block_id] = left_extent+right_extent+seed_size;
                    d_done[block_id] = true;
                    total_anchors++;
                }
            }
            else{
                d_r_starts[block_id] = ref_loc;
                d_q_starts[block_id] = query_loc;
                d_len[block_id] = seed_size;
                d_done[block_id] = false;
                total_anchors++;
            }
        }

//        if(block_id == 0 && thread_id == 0){
//            printf("%d %d %d\n", start, total_score, score[0]);
//        }

        __syncthreads();
    }
}

__global__
void find_anchors (int num_seeds, char* d_ref_seq, char* d_query_seq, uint32_t* d_index_table, uint64_t* d_pos_table, uint64_t *seed_offsets, int *d_sub_mat, int xdrop, int xdrop_threshold, uint32_t* d_r_starts, uint32_t* d_q_starts, uint32_t* d_len, bool* d_done, int ref_len, int query_len, int seed_size){
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

    __shared__ int ref_loc;
    __shared__ int query_loc;
    __shared__ int total_score;
    __shared__ int warp_size;
    int max_score;
    __shared__ bool right_xdrop_found; 
    __shared__ bool left_xdrop_found; 
    __shared__ uint32_t left_extent;
    __shared__ uint32_t right_extent;
    int thread_score;
    int current_id;
    int temp;

    __shared__ uint32_t r_starts[HIT_LIMIT];
    __shared__ uint32_t q_starts[HIT_LIMIT];
    __shared__ uint32_t hit_num[NUM_THREADS];
    __shared__ int score[NUM_THREADS];
    __shared__ uint32_t final_hits;
    __shared__ uint32_t total_hits;
    __shared__ uint32_t total_curr_hits;
    __shared__ int sub_mat[25];
    __shared__ uint32_t total_anchors;

    if(thread_id < 25){
        sub_mat[thread_id] = d_sub_mat[thread_id];
    }

    if(thread_id == 0){
        total_hits = 0;
        final_hits = 0;
        total_anchors = 0;
        total_curr_hits = 0;
        warp_size = warpSize;
    }
    int lane_id = thread_id % warp_size;
    __syncthreads();

    for (int id = id_start; id < num_seeds; id += stride) {

        current_id = id + thread_id;

        if(current_id < num_seeds){
            seed_offset = seed_offsets[current_id];

            seed = (seed_offset >> 32);
            q_start = ((seed_offset << 32) >> 32);

            // start and end from the seed block_id table
            end = d_index_table[seed];
            start = 0;
            if (seed > 0){
                start = d_index_table[seed-1];
            }
        }
        else{
            start = 0;
            end = 0;
        }

        hit_num[thread_id] = end-start;
        __syncthreads();

        // double buffer parallel prefix sum
        int k_val;
        for(int k = 1; k < block_dim; k=k*2){
            if(thread_id >= k){
                k_val = hit_num[thread_id-k];
            }
            __syncthreads();

            if(thread_id >= k){
                hit_num[thread_id] += k_val;
            }
            __syncthreads();
        }

        if(thread_id == 0){
            total_hits = hit_num[block_dim-1];
            final_hits += total_hits; 
        }
        __syncthreads();

        for(int hit_limit = 0; hit_limit < total_hits; hit_limit=hit_limit+HIT_LIMIT){
            int addr_start = (thread_id == 0) ? 0 : hit_num[thread_id-1];
            for (uint32_t p = start; p < end; p++) { 
                int index_el = addr_start+p-start;
                if ((index_el >= hit_limit) && (index_el < (hit_limit+HIT_LIMIT))) { 
                    r_starts[index_el-hit_limit] = d_pos_table[p];
                    q_starts[index_el-hit_limit] = q_start;
                }
            }
            __syncthreads();

            if(thread_id == 0){
                if(total_hits > hit_limit+HIT_LIMIT){
                    total_curr_hits = HIT_LIMIT;
                }
                else{
                    total_curr_hits = total_hits-hit_limit;
                }
            }
            __syncthreads();

            for (int id1 = 0; id1 < total_curr_hits; id1 += 1) {
                ref_loc   = r_starts[id1];
                query_loc = q_starts[id1];
                total_score = ref_loc + query_loc;

                //////////////////////////////////////////////////////////////////
                
                thread_score = 0;
                if(thread_id < seed_size){
                    thread_score = sub_mat[d_ref_seq[ref_loc+thread_id]*5+d_query_seq[query_loc+thread_id]];
                }

                for (int offset = warp_size/2; offset > 0; offset /= 2)
                    thread_score += __shfl_down_sync(0x13, thread_score, offset);

                if(thread_id == 0){
                    total_score += thread_score;
                }

                //////////////////////////////////////////////////////////////////

                thread_score = 0;
                if(ref_loc+seed_size+thread_id < ref_len && query_loc+seed_size+thread_id < query_len){
                    thread_score = sub_mat[d_ref_seq[ref_loc+seed_size+thread_id]*5+d_query_seq[query_loc+seed_size+thread_id]];
                }

                for (int offset = 1; offset < warp_size; offset *= 2){
                    temp = __shfl_up_sync(0xFFFFFFFF, thread_score, offset);

                    if(lane_id >= offset){
                        thread_score += temp;
                    }
                }
                __syncthreads();

                score[thread_id] = thread_score;
                __syncthreads();

                if(thread_id >= warp_size){
                    score[thread_id] += score[warp_size-1];
                }
                __syncthreads();

                if(thread_id == 0){
                    right_xdrop_found = false;
                    max_score = 0;
                    for(int i = 0; i < NUM_THREADS; i++){
                        if(score[i] > max_score){
                           max_score = score[i];
                        }

                        if(max_score-score[i] > xdrop && right_xdrop_found == false){
                            total_score+=max_score;
                            right_xdrop_found = true;
                            right_extent = i;
                            break;
                        }
                    }
                }
                __syncthreads();

                //////////////////////////////////////////////////////////////////

                thread_score = 0;
                if(ref_loc-thread_id-1 >= 0 && query_loc-thread_id-1 >=0){
                    thread_score = sub_mat[d_ref_seq[ref_loc-thread_id-1]*5+d_query_seq[query_loc-thread_id-1]];
                }

                for (int offset = 1; offset < warp_size; offset *= 2){
                    temp = __shfl_up_sync(0xFFFFFFFF, thread_score, offset);

                    if(lane_id >= offset){
                        thread_score += temp;
                    }
                }

                score[thread_id] = thread_score;
                __syncthreads();

                if(thread_id >= warp_size){
                    score[thread_id] += score[warp_size-1];
                }
                __syncthreads();

                if(thread_id == 0){
                    left_xdrop_found = false;
                    max_score = 0;
                    for(int i = 0; i < NUM_THREADS; i++){
                        if(score[i] > max_score){
                           max_score = score[i];
                        }

                        if(max_score-score[i] > xdrop && left_xdrop_found == false){
                            total_score+=max_score;
                            left_xdrop_found = true;
                            left_extent = i;
                            break;
                        }
                    }
                }
                __syncthreads();

                //////////////////////////////////////////////////////////////////

                if(total_anchors < OUT_LIMIT){
                    if(right_xdrop_found && left_xdrop_found){
                        if(total_score >= xdrop_threshold){
                            d_r_starts[block_id*OUT_LIMIT + total_anchors] = ref_loc - left_extent;
                            d_q_starts[block_id*OUT_LIMIT + total_anchors] = query_loc - left_extent;
                            d_len[block_id*OUT_LIMIT + total_anchors] = left_extent+right_extent+seed_size;
                            d_done[block_id*OUT_LIMIT + total_anchors] = true;
                            total_anchors++;
                        }
                    }
                    else{
                        d_r_starts[block_id*OUT_LIMIT + total_anchors] = ref_loc;
                        d_q_starts[block_id*OUT_LIMIT + total_anchors] = query_loc;
                        d_len[block_id*OUT_LIMIT + total_anchors] = seed_size;
                        d_done[block_id*OUT_LIMIT + total_anchors] = false;
                        total_anchors++;
                    }
                }
            }
            __syncthreads();
        }
    }
}

int SeedAndFilter (std::vector<uint64_t> seed_offset_vector, int n){
    int ret = 0;

    uint32_t num_seeds = seed_offset_vector.size();
    assert(num_seeds <= 13*cfg.chunk_size);
    if (num_seeds == 0) {
        return ret;
    }
    seed_size = 19;

    cudaError_t err;
    uint32_t* h_r_starts = (uint32_t*) calloc(OUT_LIMIT*NUM_BLOCKS, sizeof(uint32_t));
    uint32_t* h_q_starts = (uint32_t*) calloc(OUT_LIMIT*NUM_BLOCKS, sizeof(uint32_t));
    uint32_t* h_len      = (uint32_t*) calloc(OUT_LIMIT*NUM_BLOCKS, sizeof(uint32_t));
    bool* h_done             = (bool*) calloc(OUT_LIMIT*NUM_BLOCKS, sizeof(bool));

    for (uint32_t i = 0; i < num_seeds; i++) {
        h_seed_offsets[i] = seed_offset_vector[i];
    }

    err = cudaMemcpy(d_seed_offsets, h_seed_offsets, num_seeds*sizeof(uint64_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed!\n");
        exit(1);
    }
    
    err = cudaMemcpy(d_r_starts, h_r_starts, OUT_LIMIT*NUM_BLOCKS*sizeof(uint32_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed!\n");
        exit(1);
    }
    
    err = cudaMemcpy(d_q_starts, h_q_starts, OUT_LIMIT*NUM_BLOCKS*sizeof(uint32_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed! %d\n", err);
        exit(1);
    }
    
    err = cudaMemcpy(d_len, h_len, OUT_LIMIT*NUM_BLOCKS*sizeof(uint32_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed!\n");
        exit(1);
    }

    err = cudaMemcpy(d_done, h_done, OUT_LIMIT*NUM_BLOCKS*sizeof(bool), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed!\n");
        exit(1);
    }

//    printf("Start find_anchors %d\n", num_seeds);
//    find_anchors <<<NUM_BLOCKS, NUM_THREADS>>> (num_seeds, d_ref_seq, d_query_seq, d_index_table, d_pos_table, d_seed_offsets, d_sub_mat, cfg.xdrop, cfg.xdrop_threshold, d_r_starts, d_q_starts, d_len, d_done, ref_len, query_len, seed_size);
//    find_anchors1 <<<num_seeds, NUM_THREADS>>> (num_seeds, d_ref_seq, d_query_seq, d_index_table, d_pos_table, d_seed_offsets, d_sub_mat, cfg.xdrop, cfg.xdrop_threshold, d_r_starts, d_q_starts, d_len, d_done, ref_len, query_len, seed_size);
    find_anchors2 <<<num_seeds, NUM_THREADS>>> (num_seeds, d_ref_seq, d_query_seq, d_index_table, d_pos_table, d_seed_offsets, d_sub_mat, cfg.xdrop, cfg.xdrop_threshold, d_r_starts, d_q_starts, d_len, d_done, ref_len, query_len, seed_size);

    err = cudaMemcpy(h_r_starts, d_r_starts, OUT_LIMIT*NUM_BLOCKS*sizeof(uint32_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed!\n");
        exit(1);
    }
    
    err = cudaMemcpy(h_q_starts, d_q_starts, OUT_LIMIT*NUM_BLOCKS*sizeof(uint32_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed!\n");
        exit(1);
    }
    
    err = cudaMemcpy(h_len, d_len, OUT_LIMIT*NUM_BLOCKS*sizeof(uint32_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed!\n");
        exit(1);
    }

    err = cudaMemcpy(h_done, d_done, OUT_LIMIT*NUM_BLOCKS*sizeof(bool), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed!\n");
        exit(1);
    }

    int total_anchors = 0;
    int total_not_done = 0;

    for(int i = 0; i < OUT_LIMIT*NUM_BLOCKS; i++){
        if(h_len[i] > 0){
            if(h_done[i]){
                total_anchors++;
            }
            else{
                total_not_done++;
            }
        }
    }
//    printf("%d %d\n", total_anchors, total_not_done);

    free(h_r_starts);
    free(h_q_starts);
    free(h_len);
    free(h_done);

    gpu_lock.unlock();

    return ret;
}

size_t InitializeProcessor (int t, int f){
    size_t ret = 0;
    cudaError_t err;

    err = cudaMalloc(&d_seed_offsets, 13*cfg.chunk_size*sizeof(uint64_t)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "1 Error: cudaMalloc failed!\n");
        exit(1);
    }

    err = cudaMalloc(&d_r_starts, OUT_LIMIT*NUM_BLOCKS*sizeof(uint32_t)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "2 Error: cudaMalloc failed!\n");
        exit(1);
    }

    err = cudaMalloc(&d_q_starts, OUT_LIMIT*NUM_BLOCKS*sizeof(uint32_t)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "3 Error: cudaMalloc failed!\n");
        exit(1);
    }

    err = cudaMalloc(&d_len, OUT_LIMIT*NUM_BLOCKS*sizeof(uint32_t)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "4 Error: cudaMalloc failed!\n");
        exit(1);
    }

    err = cudaMalloc(&d_done, OUT_LIMIT*NUM_BLOCKS*sizeof(bool)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "5 Error: cudaMalloc failed!\n");
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
    
    int blockSize = 1024;
    int numBlocks = (1 << 10);
    
    compress_string <<<numBlocks, blockSize>>> (len, d_ref_seq_tmp, d_ref_seq);
    
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

    int blockSize = 1024;
    int numBlocks = (1 << 10);
    
    compress_string <<<numBlocks, blockSize>>> (len, d_query_seq_tmp, d_query_seq);

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
    cudaFree(d_r_starts);
    cudaFree(d_q_starts);
    cudaFree(d_len);
    cudaFree(d_done);
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
