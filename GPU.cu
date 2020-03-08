#include "GPU.h"
#include <thrust/scan.h>
#include <thrust/find.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>

std::mutex gpu_lock;

int err;                            
int check_status = 0;
int NUM_DEVICES;
std::vector<int> available_gpus;
std::mutex mu;
std::condition_variable cv;

uint32_t ref_len;
uint32_t query_length[BUFFER_DEPTH];

char** d_ref_seq;
char** d_query_seq;
char** d_query_rc_seq;

int **d_sub_mat;

uint32_t** d_index_table;
uint32_t** d_pos_table;

uint64_t** d_seed_offsets;

hsp** d_hsp;
hsp** d_hsp_reduced;

std::vector<thrust::device_vector<uint32_t> > d_done_vec;
std::vector<thrust::device_vector<uint32_t> > d_hit_num_vec;
uint32_t** d_done_array;
uint32_t** d_hit_num_array;

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
        else if ((ch == 'a') || (ch == 'c') || (ch == 'g') || (ch == 't')){
            dst = L_NT;
            dst_rc = L_NT;
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
        else if ((ch == 'a') || (ch == 'c') || (ch == 'g') || (ch == 't'))
            dst = L_NT;
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
void find_num_hits (int num_seeds, const uint32_t* __restrict__ d_index_table, uint64_t* seed_offsets, uint32_t* seed_hit_num){

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
void find_anchors1 (int num_seeds, const char* __restrict__  d_ref_seq, const char* __restrict__  d_query_seq, const uint32_t* __restrict__  d_index_table, const uint32_t* __restrict__ d_pos_table, uint64_t*  d_seed_offsets, int *d_sub_mat, int xdrop, int hspthresh, uint32_t* d_done, uint32_t ref_len, uint32_t query_len, int seed_size, uint32_t* seed_hit_num, int num_hits, hsp* d_hsp){

    int thread_id = threadIdx.x;
    int block_id = blockIdx.x;
    int warp_size = warpSize;
    int lane_id = thread_id%warp_size;
    int warp_id = (thread_id-lane_id)/warp_size;

    __shared__ uint32_t start, end;
    __shared__ uint32_t seed;
    __shared__ uint64_t seed_offset;

    __shared__ int ref_loc[NUM_WARPS];
    __shared__ int query_loc;
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

    __shared__ int sub_mat[NUC2];

    if(thread_id < NUC2){
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
                    thread_score = sub_mat[d_ref_seq[ref_pos]*NUC+d_query_seq[query_pos]];
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
                    thread_score = sub_mat[d_ref_seq[ref_pos]*NUC+d_query_seq[query_pos]];
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

                if(total_score[warp_id] >= hspthresh){
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

__global__
void find_anchors2 (int num_seeds, const char* __restrict__  d_ref_seq, const char* __restrict__  d_query_seq, const uint32_t* __restrict__  d_index_table, const uint32_t* __restrict__ d_pos_table, uint64_t*  d_seed_offsets, int *d_sub_mat, int xdrop, int hspthresh, uint32_t* d_done, uint32_t ref_len, uint32_t query_len, int seed_size, uint32_t* seed_hit_num, int num_hits, hsp* d_hsp){

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
    __shared__ uint32_t seed_hit_prefix;

    int thread_score;
    int max_thread_score;
    bool xdrop_done;
    int temp;
    uint32_t ref_pos;
    uint32_t query_pos;
    uint32_t pos_offset;

    __shared__ int sub_mat[NUC2];

    if(thread_id < NUC2){
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

            while(!right_xdrop_found[warp_id] && !right_edge[warp_id]){
                pos_offset = lane_id + tile[warp_id]*warp_size;
                ref_pos   = ref_loc[warp_id] + pos_offset;
                query_pos = query_loc + pos_offset;
                thread_score = 0;

                if(ref_pos < ref_len && query_pos < query_len){
                    thread_score = sub_mat[d_ref_seq[ref_pos]*NUC+d_query_seq[query_pos]];
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
                }
                else{
                    max_thread_score = prev_max_score[warp_id];
                }

                __syncwarp();
#pragma unroll
                for (int offset = 1; offset < warp_size; offset = offset << 1){
                    temp = __shfl_up_sync(0xFFFFFFFF, max_thread_score, offset);

                    if(lane_id >= offset){
                        if(temp > max_thread_score){
                            max_thread_score = temp;
                        }
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
                        right_extent[warp_id] = (tile[warp_id]+1)*warp_size;//max_pos;
                    }
                    else if(ref_pos >= ref_len || query_pos >= query_len){
                        total_score[warp_id]+=max_thread_score;
                        right_edge[warp_id] = true;
//                        right_extent[warp_id] = max_pos;
                    }
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
                pos_offset = lane_id+1+tile[warp_id]*warp_size;
                thread_score = 0;

                if(ref_loc[warp_id] >= pos_offset  && query_loc >= pos_offset){
                    ref_pos   = ref_loc[warp_id] - pos_offset;
                    query_pos = query_loc - pos_offset;
                    thread_score = sub_mat[d_ref_seq[ref_pos]*NUC+d_query_seq[query_pos]];
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
                }
                else{
                    max_thread_score = prev_max_score[warp_id];
                }
                __syncwarp();

#pragma unroll
                for (int offset = 1; offset < warp_size; offset = offset << 1){
                    temp = __shfl_up_sync(0xFFFFFFFF, max_thread_score, offset);

                    if(lane_id >= offset){
                        if(temp > max_thread_score){
                            max_thread_score = temp;
                        }
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
                        left_extent[warp_id] = (tile[warp_id]+1)*warp_size; //max_pos+1;
                    }
                    else if(ref_loc[warp_id] < pos_offset && query_loc < pos_offset){
                        total_score[warp_id]+=max_thread_score;
                        left_edge[warp_id] = true;
//                        left_extent[warp_id] = max_pos+1;
                    }
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

                if(total_score[warp_id] >= hspthresh){
                    d_hsp[dram_address].ref_start = ref_loc[warp_id] - left_extent[warp_id];
                    d_hsp[dram_address].query_start = query_loc - left_extent[warp_id];
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

std::vector<hsp> SeedAndFilter (std::vector<uint64_t> seed_offset_vector, bool rev, uint32_t buffer, uint32_t seed_size, int xdrop, int hspthresh){

    cudaError_t err;

    uint32_t num_hits;
    uint32_t num_anchors;

    uint32_t num_seeds = seed_offset_vector.size();
    assert(num_seeds <= 13*WGA_CHUNK);

uint64_t* tmp_offset = (uint64_t*) malloc(num_seeds*sizeof(uint64_t));
    for (uint32_t i = 0; i < num_seeds; i++) {
        tmp_offset[i] = seed_offset_vector[i];
    }

    hsp* h_hsp;

    int g;
    std::unique_lock<std::mutex> locker(mu);
    if (available_gpus.empty()) {
        cv.wait(locker, [](){return !available_gpus.empty();});
    }
    g = available_gpus.back();
    available_gpus.pop_back();
    locker.unlock();

    err = cudaSetDevice(g);

    err = cudaMemcpy(d_seed_offsets[g], tmp_offset, num_seeds*sizeof(uint64_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed! seed_offset\n");
        exit(1);
    }

    find_num_hits <<<MAX_BLOCKS, MAX_THREADS>>> (num_seeds, d_index_table[g], d_seed_offsets[g], d_hit_num_array[g]);

    thrust::inclusive_scan(d_hit_num_vec[g].begin(), d_hit_num_vec[g].begin() + num_seeds, d_hit_num_vec[g].begin());

    err = cudaMemcpy(&num_hits, d_hit_num_array[g]+num_seeds-1, sizeof(uint32_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed! num_hits\n");
        exit(1);
    }

//    if(rev)
//        find_anchors1 <<<num_seeds, BLOCK_SIZE>>> (num_seeds, d_ref_seq[g], d_query_rc_seq[g], d_index_table[g], d_pos_table[g], d_seed_offsets[g], d_sub_mat[g], xdrop, hspthresh, d_done[g], ref_len, query_len, seed_size, seed_hit_num_array, num_hits, d_hsp[g]);
//    else
//        find_anchors1 <<<num_seeds, BLOCK_SIZE>>> (num_seeds, d_ref_seq[g], d_query_seq[g], d_index_table[g], d_pos_table[g], d_seed_offsets[g], d_sub_mat[g], xdrop, hspthresh, d_done[g], ref_len, query_len, seed_size, seed_hit_num_array, num_hits, d_hsp[g]);

    if(rev)
        find_anchors2 <<<num_seeds, BLOCK_SIZE>>> (num_seeds, d_ref_seq[g], d_query_rc_seq[buffer*NUM_DEVICES+g], d_index_table[g], d_pos_table[g], d_seed_offsets[g], d_sub_mat[g], xdrop, hspthresh, d_done_array[g], ref_len, query_length[buffer], seed_size, d_hit_num_array[g], num_hits, d_hsp[g]);
    else
        find_anchors2 <<<num_seeds, BLOCK_SIZE>>> (num_seeds, d_ref_seq[g], d_query_seq[buffer*NUM_DEVICES+g], d_index_table[g], d_pos_table[g], d_seed_offsets[g], d_sub_mat[g], xdrop, hspthresh, d_done_array[g], ref_len, query_length[buffer], seed_size, d_hit_num_array[g], num_hits, d_hsp[g]);

    thrust::inclusive_scan(d_done_vec[g].begin(), d_done_vec[g].begin() + num_hits, d_done_vec[g].begin());

    err = cudaMemcpy(&num_anchors, d_done_array[g]+num_hits-1, sizeof(uint32_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed! num_anchors %s\n", cudaGetErrorString(err));
        exit(1);
    }

    fill_output <<<MAX_BLOCKS, MAX_THREADS>>>(d_done_array[g], d_hsp[g], d_hsp_reduced[g], num_hits);
    
    h_hsp = (hsp*) calloc(num_anchors, sizeof(hsp));

    err = cudaMemcpy(h_hsp, d_hsp_reduced[g], num_anchors*sizeof(hsp), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy failed! hsp\n");
        exit(1);
    }

    {
    std::unique_lock<std::mutex> locker(mu);
    available_gpus.push_back(g);
    locker.unlock();
    cv.notify_one();
    }

    std::vector<hsp> gpu_filter_output;
    for(int i = 0; i < num_anchors; i++){
        gpu_filter_output.push_back(h_hsp[i]);
    }
    
    free(h_hsp);
    free(tmp_offset);

    return gpu_filter_output;
}

size_t InitializeProcessor (int* sub_mat){

    size_t ret = 0;
    cudaError_t err;
    int nDevices;

    err = cudaGetDeviceCount(&nDevices);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: No GPU device found!\n");
        exit(1);
    }

    NUM_DEVICES = nDevices; 
    fprintf(stderr, "Using %d GPU(s)\n", NUM_DEVICES);

    d_seed_offsets = (uint64_t**) malloc(NUM_DEVICES*sizeof(uint64_t*));
    d_done_array = (uint32_t**) malloc(NUM_DEVICES*sizeof(uint32_t*));
    d_hit_num_array = (uint32_t**) malloc(NUM_DEVICES*sizeof(uint32_t*));
    d_hsp = (hsp**) malloc(NUM_DEVICES*sizeof(hsp*));
    d_hsp_reduced = (hsp**) malloc(NUM_DEVICES*sizeof(hsp*));
    d_sub_mat = (int**) malloc(NUM_DEVICES*sizeof(int*));
    d_done_vec.reserve(NUM_DEVICES);
    d_hit_num_vec.reserve(NUM_DEVICES);

    for(int g = 0; g < NUM_DEVICES; g++){

        cudaSetDevice(g);

        d_done_vec.emplace_back(MAX_HITS, 0);
        d_hit_num_vec.emplace_back(MAX_SEEDS, 0);

        d_done_array[g] = thrust::raw_pointer_cast(d_done_vec.at(g).data());
        d_hit_num_array[g] = thrust::raw_pointer_cast(d_hit_num_vec.at(g).data());

        err = cudaMalloc(&d_seed_offsets[g], 13*WGA_CHUNK*sizeof(uint64_t)); 
        if (err != cudaSuccess) {
            fprintf(stderr, "Error: cudaMalloc failed1!\n");
            exit(1);
        }

        err = cudaMalloc(&d_hsp[g], MAX_HITS*sizeof(hsp));
        if (err != cudaSuccess) {
            fprintf(stderr, "Error: cudaMalloc failed2!\n");
            exit(1);
        }

        err = cudaMalloc(&d_hsp_reduced[g], MAX_HITS*sizeof(hsp));
        if (err != cudaSuccess) {
            fprintf(stderr, "Error: cudaMalloc failed3!\n");
            exit(1);
        }

        err = cudaMalloc(&d_sub_mat[g], NUC2*sizeof(int)); 
        if (err != cudaSuccess) {
            fprintf(stderr, "Error: cudaMalloc failed4!\n");
            exit(1);
        }

        err = cudaMemcpy(d_sub_mat[g], sub_mat, NUC2*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) {
            fprintf(stderr, "Error: cudaMemcpy failed5!\n");
            exit(1);
        }

        available_gpus.push_back(g);
    }
    
    d_query_seq    = (char**) malloc(BUFFER_DEPTH*NUM_DEVICES*sizeof(char*));
    d_query_rc_seq = (char**) malloc(BUFFER_DEPTH*NUM_DEVICES*sizeof(char*));
    d_ref_seq = (char**) malloc(NUM_DEVICES*sizeof(char*));
    d_index_table = (uint32_t**) malloc(NUM_DEVICES*sizeof(uint32_t*));
    d_pos_table = (uint32_t**) malloc(NUM_DEVICES*sizeof(uint32_t*));

    return ret;
}

void SendRefWriteRequest (size_t start_addr, size_t len){
    cudaError_t err;
    ref_len = len;
    
    for(int g = 0; g < NUM_DEVICES; g++){

        cudaSetDevice(g);
        char* d_ref_seq_tmp;
        err = cudaMalloc(&d_ref_seq_tmp, len*sizeof(char)); 
        if (err != cudaSuccess) {
            fprintf(stderr, "Error: cudaMalloc failed6!\n");
            exit(1);
        }

        err = cudaMemcpy(d_ref_seq_tmp, g_DRAM->buffer + start_addr, len*sizeof(char), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) {
            fprintf(stderr, "Error: cudaMemcpy failed7!\n");
            exit(1);
        }

        err = cudaMalloc(&d_ref_seq[g], len*sizeof(char)); 
        if (err != cudaSuccess) {
            fprintf(stderr, "Error: cudaMalloc failed8!\n");
            exit(1);
        }

        compress_string <<<MAX_BLOCKS, MAX_THREADS>>> (len, d_ref_seq_tmp, d_ref_seq[g]);

        cudaFree(d_ref_seq_tmp);
    }
}

void SendQueryWriteRequest (size_t start_addr, size_t len, uint32_t buffer){
    cudaError_t err;
    query_length[buffer] = len;

    for(int g = 0; g < NUM_DEVICES; g++){

        cudaSetDevice(g);
        char* d_query_seq_tmp;
        err = cudaMalloc(&d_query_seq_tmp, len*sizeof(char)); 
        if (err != cudaSuccess) {
            fprintf(stderr, "Error: cudaMalloc failed9!\n");
            exit(1);
        }

        err = cudaMemcpy(d_query_seq_tmp, g_DRAM->buffer + start_addr, len*sizeof(char), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) {
            fprintf(stderr, "Error: cudaMemcpy failed!\n");
            exit(1);
        }

        err = cudaMalloc(&d_query_seq[buffer*NUM_DEVICES+g], len*sizeof(char)); 
        if (err != cudaSuccess) {
            fprintf(stderr, "Error: cudaMalloc failed10!\n");
            exit(1);
        }

        err = cudaMalloc(&d_query_rc_seq[buffer*NUM_DEVICES+g], len*sizeof(char)); 
        if (err != cudaSuccess) {
            fprintf(stderr, "Error: cudaMalloc failed11!\n");
            exit(1);
        }

        compress_string_rev_comp <<<MAX_BLOCKS, MAX_THREADS>>> (len, d_query_seq_tmp, d_query_seq[buffer*NUM_DEVICES+g], d_query_rc_seq[buffer*NUM_DEVICES+g]);

        cudaFree(d_query_seq_tmp);
    }
}

void SendSeedPosTable (uint32_t* index_table, uint32_t index_table_size, uint32_t* pos_table, uint32_t num_index){
    cudaError_t err;

    for(int g = 0; g < NUM_DEVICES; g++){

        cudaSetDevice(g);

        err = cudaMalloc(&d_index_table[g], index_table_size*sizeof(uint32_t)); 
        if (err != cudaSuccess) {
            fprintf(stderr, "Error: cudaMalloc failed!s1\n");
            exit(1);
        }

        err = cudaMemcpy(d_index_table[g], index_table, index_table_size*sizeof(uint32_t), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) {
            fprintf(stderr, "Error: cudaMemcpy failed!\n");
            exit(1);
        }

        err = cudaMalloc(&d_pos_table[g], num_index*sizeof(uint32_t)); 
        if (err != cudaSuccess) {
            fprintf(stderr, "Error: cudaMalloc failed!s2\n");
            exit(1);
        }

        err = cudaMemcpy(d_pos_table[g], pos_table, num_index*sizeof(uint32_t), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) {
            fprintf(stderr, "Error: cudaMemcpy failed!\n");
            exit(1);
        }
    }
}

void clearRef(){

    for(int g = 0; g < NUM_DEVICES; g++){

        cudaSetDevice(g);
        cudaFree(d_ref_seq[g]);
        cudaFree(d_index_table[g]);
        cudaFree(d_pos_table[g]);
    }
}

void clearQuery(uint32_t buffer){

    for(int g = 0; g < NUM_DEVICES; g++){

        cudaSetDevice(g);
        cudaFree(d_query_seq[buffer*NUM_DEVICES+g]);
        cudaFree(d_query_rc_seq[buffer*NUM_DEVICES+g]);

    }
}

void ShutdownProcessor(){

    d_done_vec.clear();
    d_hit_num_vec.clear();

    cudaDeviceReset();
}

DRAM *g_DRAM = nullptr;

InitializeProcessor_ptr g_InitializeProcessor = InitializeProcessor;
SendSeedPosTable_ptr g_SendSeedPosTable = SendSeedPosTable;
SeedAndFilter_ptr g_SeedAndFilter = SeedAndFilter;
SendRefWriteRequest_ptr g_SendRefWriteRequest = SendRefWriteRequest;
SendQueryWriteRequest_ptr g_SendQueryWriteRequest = SendQueryWriteRequest;
ShutdownProcessor_ptr g_ShutdownProcessor = ShutdownProcessor;
clearRef_ptr g_clearRef = clearRef;
clearQuery_ptr g_clearQuery = clearQuery;
