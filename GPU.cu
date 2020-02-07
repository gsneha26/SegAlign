#include "graph.h"
#include <thrust/scan.h>
#include <thrust/find.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include "tbb/scalable_allocator.h"

std::mutex gpu_lock;

const char A_NT = 0;
const char C_NT = 1;
const char G_NT = 2;
const char T_NT = 3;
const char N_NT = 4;
const char X_NT = 4;

int mat_offset[] = {0, 1, 3, 6, 10};
int err;                            
int check_status = 0;

int ref_len;
int query_len;
int seed_size;

char* d_ref_seq;
char* d_query_seq;
char* d_query_rc_seq;

int *d_sub_mat;

uint32_t* d_index_table;
uint64_t* d_pos_table;

uint64_t* h_seed_offsets;
uint64_t* d_seed_offsets;

hsp* d_hsp;
hsp* d_hsp_reduced;

thrust::device_vector<uint32_t> d_done_vec(MAX_HITS);
uint32_t* d_done = thrust::raw_pointer_cast(&d_done_vec[0]);

thrust::device_vector<int> seed_hit_num(MAX_SEEDS);
int* seed_hit_num_array = thrust::raw_pointer_cast(&seed_hit_num[0]);

__global__
void compress_string_rev_comp (uint32_t len, char* src_seq, char* dst_seq, char* dst_seq_rc){ 
    int thread_id = threadIdx.x;
    int block_dim = blockDim.x;
    int grid_dim = gridDim.x;
    int block_id = blockIdx.x;

    int stride = block_dim * grid_dim;
    uint32_t start = block_dim * block_id + thread_id;

    for (uint32_t i = start; i < len; i += stride) {
        char ch = src_seq[i];
        char dst = X_NT;
        char dst_rc = X_NT;
        if (ch == 'A'){
            dst = A_NT;
            dst_rc = T_NT;
        }
        else if (ch == 'C'){ 
            dst = C_NT;
            dst_rc = G_NT;
        }
        else if (ch == 'G'){
            dst = G_NT;
            dst_rc = C_NT;
        }
        else if (ch == 'T'){
            dst = T_NT;
            dst_rc = A_NT;
        }
        else if ((ch == 'n') || (ch == 'N')){
            dst = N_NT;
            dst_rc = N_NT;
        }
        dst_seq[i] = dst;
        dst_seq_rc[len -1 -i] = dst_rc;
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
        char dst = X_NT;
        if (ch == 'A')
            dst = A_NT;
        else if (ch == 'C')
            dst = C_NT;
        else if (ch == 'G')
            dst = G_NT;
        else if (ch == 'T')
            dst = T_NT;
        else if ((ch == 'n') || (ch == 'N'))
            dst = N_NT;
        dst_seq[i] = dst;
    }
}

__global__
void fill_output (uint32_t* d_done, hsp* d_hsp, hsp* d_hsp_reduced, int num_hits){

    int thread_id = threadIdx.x;
    int block_dim = blockDim.x;
    int grid_dim = gridDim.x;
    int block_id = blockIdx.x;

    int stride = block_dim * grid_dim;
    uint32_t start = block_dim * block_id + thread_id;
    int index = 0;

    for (uint32_t id = start; id < num_hits; id += stride) {
        index = d_done[id];

        if(id > 0){
            if(index > d_done[id-1]){
                d_hsp_reduced[index-1]    =  d_hsp[id];
            }
        }
        else{
            if(index == 1){
                d_hsp_reduced[0]    = d_hsp[0];
            }
        }
    }
}

__global__
void find_num_hits (int num_seeds, const uint32_t* __restrict__ d_index_table, uint64_t* seed_offsets, int* seed_hit_num){

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

        seed_hit_num[id] = num_seed_hit;
    }
}

__global__
void find_anchors (int num_seeds, const char* __restrict__  d_ref_seq, const char* __restrict__  d_query_seq, const uint32_t* __restrict__  d_index_table, const uint64_t* __restrict__ d_pos_table, uint64_t*  d_seed_offsets, int *d_sub_mat, int xdrop, int xdrop_threshold, uint32_t* d_done, int ref_len, int query_len, int seed_size, int* seed_hit_num, int num_hits, hsp* d_hsp){

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
    __shared__ uint32_t seed_hit_prefix;

    int thread_score;
    int max_thread_score;
    bool xdrop_done;
    int temp;
    int ref_pos;
    int query_pos;

    __shared__ int sub_mat[36];

    if(thread_id < 36){
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
        seed_hit_prefix = seed_hit_num[block_id]; 
    }

    if(warp_id == 0)
        seed_query[lane_id] = d_query_seq[query_loc+lane_id]; 
    __syncthreads();

    for (int id1 = start; id1 < end; id1 += NUM_WARPS) {
        if(id1+warp_id < end){ 
            if(lane_id == 0){ 
                ref_loc[warp_id]   = d_pos_table[id1+warp_id];
                total_score[warp_id] = 0; 
            }

            //////////////////////////////////////////////////////////////////

            thread_score = 0;
            if(lane_id < seed_size){
                thread_score = sub_mat[d_ref_seq[ref_loc[warp_id]+lane_id]*6+seed_query[lane_id]];
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
                thread_score = 0;

                if(ref_pos < ref_len && query_pos < query_len){
                    thread_score = sub_mat[d_ref_seq[ref_pos]*6+d_query_seq[query_pos]];
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
                        right_extent[warp_id] = (tile[warp_id]+1)*warp_size-1;
                        if(ref_pos >= ref_len || query_pos >= query_len)
                            right_extent[warp_id]-=max(ref_pos-ref_len+1, query_pos-query_len+1);
                    }
                    else if(ref_pos >= ref_len || query_pos >= query_len)
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
                thread_score = 0;

                if(ref_pos >= 0  && query_pos >= 0){
                    thread_score = sub_mat[d_ref_seq[ref_pos]*6+d_query_seq[query_pos]];
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
                        left_extent[warp_id] = (tile[warp_id]+1)*warp_size-1;
                        if(ref_pos < 0 || query_pos < 0)
                            right_extent[warp_id]+=min(ref_pos, query_pos);
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
                int dram_address = seed_hit_prefix -id1 - warp_id+start-1;

                if(total_score[warp_id] >= xdrop_threshold){
                    d_hsp[dram_address].ref_start = ref_loc[warp_id] - left_extent[warp_id];
                    d_hsp[dram_address].query_start = query_loc - left_extent[warp_id];
                    d_hsp[dram_address].len = left_extent[warp_id]+right_extent[warp_id]+seed_size;
                    d_hsp[dram_address].score = total_score[warp_id];
                    d_done[dram_address] = 1;

                }
                else{
                    d_hsp[dram_address].ref_start = 0;
                    d_hsp[dram_address].query_start = 0;
                    d_hsp[dram_address].len = 0;
                    d_hsp[dram_address].score = 0;
                    d_done[dram_address] = 0;
                }
            }
            __syncwarp();
        }
    }
}

__global__
void find_anchors1 (int num_seeds, const char* __restrict__  d_ref_seq, const char* __restrict__  d_query_seq, const uint32_t* __restrict__  d_index_table, const uint64_t* __restrict__ d_pos_table, uint64_t*  d_seed_offsets, int *d_sub_mat, int xdrop, int xdrop_threshold, uint32_t* d_done, int ref_len, int query_len, int seed_size, int* seed_hit_num, int num_hits, hsp* d_hsp){

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
    __shared__ int prev_max_pos[NUM_WARPS];
    __shared__ bool right_edge[NUM_WARPS]; 
    __shared__ bool left_edge[NUM_WARPS]; 
    __shared__ bool right_xdrop_found[NUM_WARPS]; 
    __shared__ bool left_xdrop_found[NUM_WARPS]; 
    __shared__ uint32_t left_extent[NUM_WARPS];
    __shared__ uint32_t right_extent[NUM_WARPS];
    __shared__ uint32_t tile[NUM_WARPS];
    __shared__ uint32_t seed_hit_prefix;

    int thread_score;
    int max_thread_score;
    bool xdrop_done;
    int temp;
    int temp_pos;
    int ref_pos;
    int query_pos;
    int max_pos;

    __shared__ int sub_mat[36];

    if(thread_id < 36){
        sub_mat[thread_id] = d_sub_mat[thread_id];
    }

    if(thread_id == 0){
        seed_offset = d_seed_offsets[block_id];
        seed = (seed_offset >> 32);
        query_loc = ((seed_offset << 32) >> 32) + seed_size;

        // start and end from the seed block_id table
        end = d_index_table[seed];
        start = 0;
        if (seed > 0){
            start = d_index_table[seed-1];
        }
        seed_hit_prefix = seed_hit_num[block_id]; 
    }
    __syncthreads();

    for (int id1 = start; id1 < end; id1 += NUM_WARPS) {
        if(id1+warp_id < end){ 
            if(lane_id == 0){ 
                ref_loc[warp_id]   = d_pos_table[id1+warp_id] + seed_size;
                total_score[warp_id] = 0; 
            }

            //////////////////////////////////////////////////////////////////

            tile[warp_id] = 0;
            right_xdrop_found[warp_id] = false;
            right_edge[warp_id] = false;
            prev_score[warp_id] = 0;
            prev_max_score[warp_id] = 0;
            right_extent[warp_id] = 0;
            prev_max_pos[warp_id] = 0;
            max_pos = 0;

            while(!right_xdrop_found[warp_id] && !right_edge[warp_id]){
                ref_pos   = ref_loc[warp_id]  + lane_id + tile[warp_id]*warp_size;
                query_pos = query_loc + lane_id + tile[warp_id]*warp_size;
                thread_score = 0;

                if(ref_pos < ref_len && query_pos < query_len){
                    thread_score = sub_mat[d_ref_seq[ref_pos]*6+d_query_seq[query_pos]];
                }

#pragma unroll
                for (int offset = 1; offset < warp_size; offset = offset << 1){
                    temp = __shfl_up_sync(0xFFFFFFFF, thread_score, offset);

                    if(lane_id >= offset){
                        thread_score += temp;
                    }
                }

                thread_score += prev_score[warp_id];
                if(thread_score > prev_max_score[warp_id]){
                     max_thread_score = thread_score;
                     max_pos = tile[warp_id]*warp_size + lane_id; 
                }
                else{
                    max_thread_score = prev_max_score[warp_id];
                    max_pos = prev_max_pos[warp_id]; 
                }

                __syncwarp();
#pragma unroll
                for (int offset = 1; offset < warp_size; offset = offset << 1){
                    temp = __shfl_up_sync(0xFFFFFFFF, max_thread_score, offset);
                    temp_pos = __shfl_up_sync(0xFFFFFFFF, max_pos, offset);

                    if(lane_id >= offset){
                        if(temp > max_thread_score){
                            max_thread_score = temp;
                            max_pos = temp_pos;
                        }
                    }
                }

//                if(ref_loc[warp_id] == 10317345 && query_loc == 6692)    
//                    printf("Right %d %d %d %d\n", thread_id, thread_score, max_thread_score, max_pos);

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
                        right_extent[warp_id] = max_pos;
                    }
                    else if(ref_pos >= ref_len || query_pos >= query_len){
                        total_score[warp_id]+=max_thread_score;
                        right_edge[warp_id] = true;
                        right_extent[warp_id] = max_pos;
                    }
                    else{
                        prev_score[warp_id] = thread_score;
                        prev_max_score[warp_id] = max_thread_score;
                        prev_max_pos[warp_id] = max_pos;
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
            prev_max_pos[warp_id] = 0;
            max_pos = 0;

            while(!left_xdrop_found[warp_id] && !left_edge[warp_id]){

                ref_pos   = ref_loc[warp_id] - lane_id - 1 - tile[warp_id]*warp_size;
                query_pos = query_loc - lane_id - 1 - tile[warp_id]*warp_size;
                thread_score = 0;

                if(ref_pos >= 0  && query_pos >= 0){
                    thread_score = sub_mat[d_ref_seq[ref_pos]*6+d_query_seq[query_pos]];
                }

#pragma unroll
                for (int offset = 1; offset < warp_size; offset = offset << 1){
                    temp = __shfl_up_sync(0xFFFFFFFF, thread_score, offset);

                    if(lane_id >= offset){
                        thread_score += temp;
                    }
                }

                thread_score += prev_score[warp_id];
                if(thread_score > prev_max_score[warp_id]){
                     max_thread_score = thread_score;
                     max_pos = tile[warp_id]*warp_size + lane_id; 
                }
                else{
                    max_thread_score = prev_max_score[warp_id];
                    max_pos = prev_max_pos[warp_id]; 
                }
                __syncwarp();

#pragma unroll
                for (int offset = 1; offset < warp_size; offset = offset << 1){
                    temp = __shfl_up_sync(0xFFFFFFFF, max_thread_score, offset);
                    temp_pos = __shfl_up_sync(0xFFFFFFFF, max_pos, offset);

                    if(lane_id >= offset){
                        if(temp > max_thread_score){
                            max_thread_score = temp;
                            max_pos = temp_pos;
                        }
                    }
                }

                xdrop_done = ((max_thread_score-thread_score) > xdrop);
                __syncwarp();

#pragma unroll
                for (int offset = 1; offset < warp_size; offset = offset << 1){
                    xdrop_done |= __shfl_up_sync(0xFFFFFFFF, xdrop_done, offset);
                }

//                if(ref_loc[warp_id] == 10317345 && query_loc == 6692)    
//                    printf("Left %d %d %d %d\n", thread_id, thread_score, max_thread_score, max_pos);

                if(lane_id == warp_size-1){
                    if(xdrop_done){
                        total_score[warp_id]+=max_thread_score;
                        left_xdrop_found[warp_id] = true;
                        left_extent[warp_id] = max_pos+1;
                    }
                    else if(ref_pos < 0 || query_pos < 0){
                        total_score[warp_id]+=max_thread_score;
                        left_edge[warp_id] = true;
                        left_extent[warp_id] = max_pos+1;
                    }
                    else{
                        prev_score[warp_id] = thread_score;
                        prev_max_score[warp_id] = max_thread_score;
                        prev_max_pos[warp_id] = max_pos;
                    }
                }
                __syncwarp();

                tile[warp_id]++;
            }

            //////////////////////////////////////////////////////////////////

            if(lane_id == 0){

                int dram_address = seed_hit_prefix -id1 - warp_id+start-1;

                if(total_score[warp_id] >= xdrop_threshold){
                    d_hsp[dram_address].ref_start = ref_loc[warp_id] - left_extent[warp_id];
                    d_hsp[dram_address].query_start = query_loc - left_extent[warp_id];
//                    d_hsp[dram_address].ref_start = ref_loc[warp_id];
//                    d_hsp[dram_address].query_start = query_loc;
                    d_hsp[dram_address].len = left_extent[warp_id]+right_extent[warp_id];
                    d_hsp[dram_address].score = total_score[warp_id];
                    d_done[dram_address] = 1;
                }
                else{
                    d_hsp[dram_address].ref_start = 0;
                    d_hsp[dram_address].query_start = 0;
                    d_hsp[dram_address].len = 0;
                    d_hsp[dram_address].score = 0;
                    d_done[dram_address] = 0;
                }
            }
            __syncwarp();
        }
    }
}

std::vector<hsp> SeedAndFilter (std::vector<uint64_t> seed_offset_vector, bool rev){

    cudaError_t err;
    seed_size = 19;

    uint32_t num_hits;
    uint32_t num_anchors;

    uint32_t num_seeds = seed_offset_vector.size();
    assert(num_seeds <= 13*cfg.chunk_size);

    hsp* h_hsp;

    gpu_lock.lock();

    for (uint32_t i = 0; i < num_seeds; i++) {
        h_seed_offsets[i] = seed_offset_vector[i];
    }

    err = cudaMemcpy(d_seed_offsets, h_seed_offsets, num_seeds*sizeof(uint64_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed!\n");
        exit(1);
    }

    find_num_hits <<<MAX_BLOCKS, MAX_THREADS>>> (num_seeds, d_index_table, d_seed_offsets, seed_hit_num_array);

    thrust::inclusive_scan(seed_hit_num.begin(), seed_hit_num.begin() + num_seeds, seed_hit_num.begin());

    err = cudaMemcpy(&num_hits, (seed_hit_num_array+num_seeds-1), sizeof(uint32_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed!\n");
        exit(1);
    }

//    if(rev)
//        find_anchors <<<num_seeds,BLOCK_SIZE>>> (num_seeds, d_ref_seq, d_query_rc_seq, d_index_table, d_pos_table, d_seed_offsets, d_sub_mat, cfg.xdrop, cfg.xdrop_threshold, d_done, ref_len, query_len, seed_size, seed_hit_num_array, num_hits, d_hsp);
//    else
//        find_anchors <<<num_seeds,BLOCK_SIZE>>> (num_seeds, d_ref_seq, d_query_seq, d_index_table, d_pos_table, d_seed_offsets, d_sub_mat, cfg.xdrop, cfg.xdrop_threshold, d_done, ref_len, query_len, seed_size, seed_hit_num_array, num_hits, d_hsp);

    if(rev)
        find_anchors1 <<<num_seeds,BLOCK_SIZE>>> (num_seeds, d_ref_seq, d_query_rc_seq, d_index_table, d_pos_table, d_seed_offsets, d_sub_mat, cfg.xdrop, cfg.xdrop_threshold, d_done, ref_len, query_len, seed_size, seed_hit_num_array, num_hits, d_hsp);
    else
        find_anchors1 <<<num_seeds,BLOCK_SIZE>>> (num_seeds, d_ref_seq, d_query_seq, d_index_table, d_pos_table, d_seed_offsets, d_sub_mat, cfg.xdrop, cfg.xdrop_threshold, d_done, ref_len, query_len, seed_size, seed_hit_num_array, num_hits, d_hsp);

    thrust::inclusive_scan(d_done_vec.begin(), d_done_vec.begin() + num_hits, d_done_vec.begin());

    err = cudaMemcpy(&num_anchors, (d_done+num_hits-1), sizeof(uint32_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed!\n");
        exit(1);
    }

    fill_output <<<MAX_BLOCKS, MAX_THREADS>>>(d_done, d_hsp, d_hsp_reduced, num_hits);
    
    h_hsp = (hsp*) calloc(num_anchors, sizeof(hsp));

    err = cudaMemcpy(h_hsp, d_hsp_reduced, num_anchors*sizeof(hsp), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed!\n");
        exit(1);
    }

    gpu_lock.unlock();

    std::vector<hsp> gpu_filter_output;
    for(int i = 0; i < num_anchors; i++){
        gpu_filter_output.push_back(h_hsp[i]);
//        if(!rev)
//            printf("%d %d %d\n", h_hsp[i].ref_start, h_hsp[i].query_start, h_hsp[i].score);
    }
    
//    printf("%d %d %d\n", num_seeds, num_hits, num_anchors);

    free(h_hsp);

    return gpu_filter_output;
}

size_t InitializeProcessor (int t, int f){
    size_t ret = 0;
    cudaError_t err;

    err = cudaMalloc(&d_seed_offsets, 13*cfg.chunk_size*sizeof(uint64_t)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMalloc failed1!\n");
        exit(1);
    }

    err = cudaMalloc(&d_hsp, MAX_HITS*sizeof(hsp));
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMalloc failed2!\n");
        exit(1);
    }

    err = cudaMalloc(&d_hsp_reduced, MAX_HITS*sizeof(hsp));
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMalloc failed3!\n");
        exit(1);
    }

    int *sub_mat;
    sub_mat = (int *)malloc(36 * sizeof(int)); 

    sub_mat[28] = cfg.gact_sub_mat[10];
    sub_mat[29] = -1000;//cfg.gact_sub_mat[10];
    sub_mat[34] = -1000;//cfg.gact_sub_mat[10];
    sub_mat[35] = -1000;//cfg.gact_sub_mat[10];
    for(int i = 0; i < 4; i++){
        sub_mat[i*6+4] = cfg.gact_sub_mat[10];
        sub_mat[4*6+i] = cfg.gact_sub_mat[10];
        sub_mat[i*6+5] = -1000;//cfg.gact_sub_mat[10];
        sub_mat[5*6+i] = -1000;//cfg.gact_sub_mat[10];
        sub_mat[i*6+i] = cfg.gact_sub_mat[i*4 + i - mat_offset[i]];
    }

    for(int i = 0; i < 4; i++){
        for(int j = i+1; j < 4; j++){
            sub_mat[i*6+j] = cfg.gact_sub_mat[i*4 + j - mat_offset[i]];
            sub_mat[j*6+i] = sub_mat[i*6+j];
        }
    }

    err = cudaMalloc(&d_sub_mat, 36*sizeof(int)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMalloc failed!\n");
        exit(1);
    }

    err = cudaMemcpy(d_sub_mat, sub_mat, 36*sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed!\n");
        exit(1);
    }

    free(sub_mat);

    return ret;
}

void SendRefWriteRequest (size_t start_addr, size_t len){
    cudaError_t err;
    ref_len = len;
    printf("Ref len: %lu\n", len);
    
    char* d_ref_seq_tmp;
    err = cudaMalloc(&d_ref_seq_tmp, len*sizeof(char)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMalloc failed!\n");
        exit(1);
    }

    err = cudaMemcpy(d_ref_seq_tmp, g_DRAM->buffer + start_addr, len*sizeof(char), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed!\n");
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
    printf("Query len: %lu\n", len);
    
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

    compress_string_rev_comp <<<MAX_BLOCKS, MAX_THREADS>>> (len, d_query_seq_tmp, d_query_seq, d_query_rc_seq);

    cudaFree(d_query_seq_tmp);
}

void SendSeedPosTable (uint32_t* index_table, uint32_t index_table_size, uint64_t* pos_table, uint32_t num_index){
    cudaError_t err;

    h_seed_offsets = (uint64_t*) malloc(13*cfg.chunk_size*sizeof(uint64_t));

    err = cudaMalloc(&d_index_table, index_table_size*sizeof(uint32_t)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMalloc failed!s1\n");
        exit(1);
    }
    
    err = cudaMemcpy(d_index_table, index_table, index_table_size*sizeof(uint32_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed!\n");
        exit(1);
    }

    err = cudaMalloc(&d_pos_table, num_index*sizeof(uint64_t)); 
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMalloc failed!s2\n");
        exit(1);
    }

    err = cudaMemcpy(d_pos_table, pos_table, num_index*sizeof(uint64_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed!\n");
        exit(1);
    }
}

void ShutdownProcessor(){
    cudaFree(d_ref_seq);
    cudaFree(d_query_seq);
    cudaFree(d_query_rc_seq);

    cudaFree(d_sub_mat);

    cudaFree(d_index_table);
    cudaFree(d_pos_table);

    free(h_seed_offsets);
    cudaFree(d_seed_offsets);

    cudaFree(d_hsp);
    cudaFree(d_hsp_reduced);
}

DRAM *g_DRAM = nullptr;

InitializeProcessor_ptr g_InitializeProcessor = InitializeProcessor;
SendSeedPosTable_ptr g_SendSeedPosTable = SendSeedPosTable;
SeedAndFilter_ptr g_SeedAndFilter = SeedAndFilter;
SendRefWriteRequest_ptr g_SendRefWriteRequest = SendRefWriteRequest;
SendQueryWriteRequest_ptr g_SendQueryWriteRequest = SendQueryWriteRequest;       
ShutdownProcessor_ptr g_ShutdownProcessor = ShutdownProcessor;
