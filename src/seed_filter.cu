#include <thrust/binary_search.h>
#include <unistd.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/scan.h>
#include <thrust/unique.h>
#include "cuda_utils.h"
#include "parameters.h"
#include "seed_filter.h"
#include "seed_filter_interface.h"
#include "store.h"
#include "store_gpu.h"

// Each ScoredSegmentPair is 16B
// With 64MB for the HSPs array per 1GB GPU memory
// With higher GPU memory, the size just linearly increases

#define MAX_HITS_PER_GB 4194304

int MAX_SEEDS;
int MAX_HITS;

int32_t seed_size;
int32_t hspthresh;

uint32_t** d_hit_num_array;
std::vector<thrust::device_vector<uint32_t> > d_hit_num_vec;

std::vector<device_buffer<uint64_t>> d_seed_offsets;
std::vector<device_buffer<ScoredSegmentPair>> d_hsp;
std::vector<device_buffer<SeedPair>> d_hit;
std::vector<DefaultDeviceAllocator> allocator_;
std::vector<CudaStream> stream_;
std::vector<std::unique_ptr<Extender>> ungapped_extender_;

std::vector<device_buffer<int8_t>> d_query_seq;
std::vector<device_buffer<int8_t>> d_query_rc_seq;
std::vector<device_buffer<int8_t>> d_target_seq;
int32_t query_length[BUFFER_DEPTH];

__global__
void compress_string_rev_comp (uint32_t len, 
                               char* src_seq, 
                               int8_t* dst_seq, 
                               int8_t* dst_seq_rc){ 
    int thread_id = threadIdx.x;
    int block_dim = blockDim.x;
    int grid_dim = gridDim.x;
    int block_id = blockIdx.x;

    int stride = block_dim * grid_dim;
    uint32_t start = block_dim * block_id + thread_id;

    for (uint32_t i = start; i < len; i += stride) {
        char ch = src_seq[i];
        int8_t dst = X_ANT;
        int8_t dst_rc = X_ANT;
        if (ch == 'A'){
            dst = A_ANT;
            dst_rc = T_ANT;
        }
        else if (ch == 'C'){ 
            dst = C_ANT;
            dst_rc = G_ANT;
        }
        else if (ch == 'G'){
            dst = G_ANT;
            dst_rc = C_ANT;
        }
        else if (ch == 'T'){
            dst = T_ANT;
            dst_rc = A_ANT;
        }
        else if ((ch == 'a') || (ch == 'c') || (ch == 'g') || (ch == 't')){
            dst = L_ANT;
            dst_rc = L_ANT;
        }
        else if ((ch == 'n') || (ch == 'N')){
            dst = N_ANT;
            dst_rc = N_ANT;
        }
        else if (ch == '&'){
            dst = E_ANT;
            dst_rc = E_ANT;
        }
        dst_seq[i] = dst;
        dst_seq_rc[len -1 -i] = dst_rc;
    }
}

__global__
void find_num_hits (int num_seeds, 
                    const uint32_t* __restrict__ d_index_table, 
                    uint64_t* seed_offsets, 
                    uint32_t* seed_hit_num){

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
void find_hits (const uint32_t* __restrict__  d_index_table, 
                const uint32_t* __restrict__ d_pos_table, 
                uint64_t*  d_seed_offsets, 
                uint32_t seed_size, 
                uint32_t* seed_hit_num, 
                int num_hits, 
                SeedPair* d_hit, 
                uint32_t start_seed_index, 
                uint32_t start_hit_index){

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
    __shared__ uint32_t seed_hit_prefix;

    if(thread_id == 0){
        seed_offset = d_seed_offsets[block_id+start_seed_index];
        seed = (seed_offset >> 32);
        query_loc = ((seed_offset << 32) >> 32) + seed_size;

        // start and end from the seed block_id table
        end = d_index_table[seed];
        start = 0;
        if (seed > 0){
            start = d_index_table[seed-1];
        }
        seed_hit_prefix = seed_hit_num[block_id+start_seed_index]; 
    }
    __syncthreads();


    for (int id1 = start; id1 < end; id1 += NUM_WARPS) {
        if(id1+warp_id < end){ 
            if(lane_id == 0){ 
                ref_loc[warp_id]   = d_pos_table[id1+warp_id] + seed_size;
                int dram_address = seed_hit_prefix -id1 - warp_id+start-1-start_hit_index;

                d_hit[dram_address].target_position_in_read = ref_loc[warp_id];
                d_hit[dram_address].query_position_in_read = query_loc; 
            }
        }
    }
}

std::vector<ScoredSegmentPair> SeedAndFilter (std::vector<uint64_t> seed_offset_vector, 
                                              bool rev, 
                                              uint32_t buffer){
    uint32_t num_hits = 0;
    uint32_t total_anchors = 0;

    uint32_t num_seeds = seed_offset_vector.size();
    assert(num_seeds <= MAX_SEEDS);

    uint64_t* tmp_offset = (uint64_t*) malloc(num_seeds*sizeof(uint64_t));
    for (uint32_t i = 0; i < num_seeds; i++) {
        tmp_offset[i] = seed_offset_vector[i];
    }

    int g;
    std::unique_lock<std::mutex> locker(mu);
    if (available_gpus.empty()) {
        cv.wait(locker, [](){return !available_gpus.empty();});
    }
    g = available_gpus.back();
    available_gpus.pop_back();
    locker.unlock();

    check_cuda_setDevice(g, "SeedAndFilter");

    check_cuda_memcpy(d_seed_offsets[g].data(), 
                      (void*)tmp_offset, 
                      num_seeds*sizeof(uint64_t), 
                      cudaMemcpyHostToDevice, 
                      "seed_offsets");

    find_num_hits <<<MAX_BLOCKS, MAX_THREADS>>> (num_seeds, 
                                                 d_index_table[g].data(), 
                                                 d_seed_offsets[g].data(), 
                                                 d_hit_num_array[g]);

    thrust::inclusive_scan(d_hit_num_vec[g].begin(),
                           d_hit_num_vec[g].begin() + num_seeds, 
                           d_hit_num_vec[g].begin());

    check_cuda_memcpy((void*)&num_hits, 
                      (void*)(d_hit_num_array[g]+num_seeds-1), 
                      sizeof(uint32_t), 
                      cudaMemcpyDeviceToHost, 
                      "num_hits");
    
    int num_iter = num_hits/MAX_HITS+1;
    uint32_t iter_hit_limit = MAX_HITS;
    thrust::device_vector<uint32_t> limit_pos (num_iter); 

    for(int i = 0; i < num_iter-1; i++){
        thrust::device_vector<uint32_t>::iterator result_end = thrust::lower_bound(d_hit_num_vec[g].begin(), 
                                                                                   d_hit_num_vec[g].begin()+num_seeds, 
                                                                                   iter_hit_limit);
        uint32_t pos = thrust::distance(d_hit_num_vec[g].begin(), result_end)-1;
        iter_hit_limit = d_hit_num_vec[g][pos]+MAX_HITS;
        limit_pos[i] = pos;
    }

    limit_pos[num_iter-1] = num_seeds-1;

    ScoredSegmentPair** h_hsp = (ScoredSegmentPair**) malloc(num_iter*sizeof(ScoredSegmentPair*));
    int32_t* num_anchors = (int32_t*) calloc(num_iter, sizeof(int32_t));

    uint32_t start_seed_index = 0;
    uint32_t start_hit_val = 0;
    uint32_t iter_num_seeds, iter_num_hits;

    if(num_hits > 0){
        
        for(int i = 0; i < num_iter; i++){
            iter_num_seeds = limit_pos[i] + 1 - start_seed_index;
            iter_num_hits  = d_hit_num_vec[g][limit_pos[i]] - start_hit_val;

            find_hits <<<iter_num_seeds, BLOCK_SIZE>>> (d_index_table[g].data(), 
                                                        d_pos_table[g].data(), 
                                                        d_seed_offsets[g].data(), 
                                                        seed_size, 
                                                        d_hit_num_array[g], 
                                                        iter_num_hits, 
                                                        d_hit[g].data(), 
                                                        start_seed_index, 
                                                        start_hit_val);

            device_buffer<int32_t> d_num_hsp(1, allocator_[g], stream_[g].get());

            if(rev){
                ungapped_extender_[g]->extend_async(d_ref_seq[g].data(),
                                                    ref_len,
                                                    d_query_rc_seq[buffer*NUM_DEVICES+g].data(),
                                                    query_length[buffer],
                                                    hspthresh,
                                                    d_hit[g].data(),
                                                    iter_num_hits,
                                                    d_hsp[g].data(),
                                                    d_num_hsp.data());
            }
            else{
                ungapped_extender_[g]->extend_async(d_ref_seq[g].data(),
                                                    ref_len,
                                                    d_query_seq[buffer*NUM_DEVICES+g].data(),
                                                    query_length[buffer],
                                                    hspthresh,
                                                    d_hit[g].data(),
                                                    iter_num_hits,
                                                    d_hsp[g].data(),
                                                    d_num_hsp.data());
            }

            num_anchors[i] = cudautils::get_value_from_device(d_num_hsp.data(), stream_[g].get());

            if(num_anchors[i] > 0){
                total_anchors += num_anchors[i];

                h_hsp[i] = (ScoredSegmentPair*) calloc(num_anchors[i], sizeof(ScoredSegmentPair));

                check_cuda_memcpy((void*)h_hsp[i], d_hsp[g].data(), num_anchors[i]*sizeof(ScoredSegmentPair), cudaMemcpyDeviceToHost, "hsp_output");
            }

            start_seed_index = limit_pos[i] + 1;
            start_hit_val = d_hit_num_vec[g][limit_pos[i]];
        }
    }

    limit_pos.clear();

    {
        std::unique_lock<std::mutex> locker(mu);
        available_gpus.push_back(g);
        locker.unlock();
        cv.notify_one();
    }
    std::vector<ScoredSegmentPair> gpu_filter_output;

    ScoredSegmentPair first_el;
    first_el.length = total_anchors;
    first_el.score = num_hits;
    gpu_filter_output.push_back(first_el);

    if(total_anchors > 0){
        for(int it = 0; it < num_iter; it++){

            for(int i = 0; i < num_anchors[it]; i++){
                gpu_filter_output.push_back(h_hsp[it][i]);
            }

            if(num_anchors[it] > 0){
                free(h_hsp[it]);
            }
        }
    }
    
    free(h_hsp);
    free(num_anchors);
    free(tmp_offset);

    return gpu_filter_output;
}

void InitializeProcessor (bool transition, 
                          uint32_t WGA_CHUNK, 
                          uint32_t input_seed_size, 
                          int* sub_mat, 
                          int input_xdrop, 
                          int input_hspthresh, 
                          bool input_noentropy){

    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, 0);
    float global_mem_gb = static_cast<float>(deviceProp.totalGlobalMem / 1073741824.0f);

    if(transition)
        MAX_SEEDS = 13*WGA_CHUNK;
    else
        MAX_SEEDS = WGA_CHUNK;

    MAX_HITS = MAX_HITS_PER_GB*global_mem_gb;

    seed_size = input_seed_size;
    hspthresh = input_hspthresh;

    d_hit_num_array = (uint32_t**) malloc(NUM_DEVICES*sizeof(uint32_t*));
    d_hit_num_vec.reserve(NUM_DEVICES);

    for(int g = 0; g < NUM_DEVICES; g++){

        check_cuda_setDevice(g, "InitializeProcessor");

        d_hit_num_vec.emplace_back(MAX_SEEDS, 0);
        d_hit_num_array[g] = thrust::raw_pointer_cast(d_hit_num_vec.at(g).data());

        stream_.push_back(make_cuda_stream());
        const int64_t max_gpu_memory     = cudautils::find_largest_contiguous_device_memory_section()/2;
        allocator_.push_back(create_default_device_allocator(max_gpu_memory, stream_[g].get()));

        d_seed_offsets.push_back(device_buffer<uint64_t>(MAX_SEEDS, allocator_[g], stream_[g].get()));
        d_hit.push_back(device_buffer<SeedPair>(MAX_HITS, allocator_[g], stream_[g].get()));
        d_hsp.push_back(device_buffer<ScoredSegmentPair>(MAX_HITS, allocator_[g], stream_[g].get()));
        ungapped_extender_.push_back(create_extender(sub_mat, NUC2, input_xdrop, input_noentropy, stream_[g].get(), g, allocator_[g]));

        available_gpus.push_back(g);
    }
}

void SendQueryWriteRequest (size_t start_addr, 
                            uint32_t len, 
                            uint32_t buffer){

    query_length[buffer] = len;

    for(int g = 0; g < NUM_DEVICES; g++){

        check_cuda_setDevice(g, "SendQueryWriteRequest");

        device_buffer<char> d_query_seq_tmp(len, 
                                            allocator_[g], 
                                            stream_[g].get());

        check_cuda_memcpy(d_query_seq_tmp.data(), 
                          (void*)(query_DRAM->buffer + start_addr), 
                          len*sizeof(int8_t), 
                          cudaMemcpyHostToDevice, 
                          "query_seq");

        d_query_seq.push_back(device_buffer<int8_t>(len, allocator_[g], stream_[g].get()));
        d_query_rc_seq.push_back(device_buffer<int8_t>(len, allocator_[g], stream_[g].get()));

        compress_string_rev_comp <<<MAX_BLOCKS, MAX_THREADS>>> (len, 
                                                                d_query_seq_tmp.data(), 
                                                                d_query_seq[buffer*NUM_DEVICES+g].data(), 
                                                                d_query_rc_seq[buffer*NUM_DEVICES+g].data());
    }
}

void ClearQuery(uint32_t buffer){

        d_query_seq.clear();
        d_query_rc_seq.clear();
        d_hit.clear();
        d_hsp.clear();
}

void ShutdownProcessor(){

    d_seed_offsets.clear();
    d_hit.clear();
    d_hit_num_vec.clear();
    d_hsp.clear();

    d_index_table.clear();
    d_pos_table.clear();
    d_ref_seq.clear();
    d_query_seq.clear();
    d_query_rc_seq.clear();

    ungapped_extender_.clear();
    allocator_.clear();
    stream_.clear();
    cudaDeviceReset();
}

InitializeProcessor_ptr g_InitializeProcessor = InitializeProcessor;
SendQueryWriteRequest_ptr g_SendQueryWriteRequest = SendQueryWriteRequest;
SeedAndFilter_ptr g_SeedAndFilter = SeedAndFilter;
ClearQuery_ptr g_ClearQuery = ClearQuery;
ShutdownProcessor_ptr g_ShutdownProcessor = ShutdownProcessor;
