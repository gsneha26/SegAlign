#include "seed_filter_interface.h"
#include "cuda_utils.h"
#include "tbb/parallel_sort.h"
#include "ntcoding.h"
#include "store_gpu.h"
#include <thrust/scan.h>
#include <thrust/execution_policy.h>

void InclusivePrefixScan (uint32_t* data, uint32_t len) {
    int g;
    
    {
        std::unique_lock<std::mutex> locker(mu);
        if (available_gpus.empty()) {
            cv.wait(locker, [](){return !available_gpus.empty();});
        }
        g = available_gpus.back();
        available_gpus.pop_back();
        locker.unlock();

        check_cuda_setDevice(g, "InclusivePrefixScan");
    }

    thrust::inclusive_scan(thrust::host, data, data + len, data); 

    {
        std::unique_lock<std::mutex> locker(mu);
        available_gpus.push_back(g);
        locker.unlock();
        cv.notify_one();
    }
}

void SendSeedPosTable (uint32_t* index_table, uint32_t index_table_size, uint32_t* pos_table, uint32_t num_index){

    for(int g = 0; g < NUM_DEVICES; g++){

        check_cuda_setDevice(g, "SendSeedPosTable");

        check_cuda_malloc((void**)&d_index_table[g], index_table_size*sizeof(uint32_t), "index_table"); 

        check_cuda_memcpy((void*)d_index_table[g], (void*)index_table, index_table_size*sizeof(uint32_t), cudaMemcpyHostToDevice, "index_table");

        check_cuda_malloc((void**)&d_pos_table[g], num_index*sizeof(uint32_t), "pos_table"); 

        check_cuda_memcpy((void*)d_pos_table[g], (void*)pos_table, num_index*sizeof(uint32_t), cudaMemcpyHostToDevice, "pos_table");
    }
}

void GenerateSeedPosTable(char* ref_str, size_t start_addr, uint32_t ref_length, uint32_t step, int shape_size, int kmer_size) {

    assert(kmer_size <= 15);
    assert(kmer_size > 3); 

    uint32_t *index_table;
    uint32_t *pos_table;
    uint32_t index_table_size;

    uint32_t offset = (shape_size+1)%step;
    uint32_t start_offset = step - offset;
    
    index_table_size = ((uint32_t)1 << 2*kmer_size) + 1;
    index_table = (uint32_t*) calloc(index_table_size, sizeof(uint32_t));

    uint32_t num_steps = (ref_length - shape_size + offset) / step;

    uint32_t* tmp_index_arr = (uint32_t*) malloc(num_steps * sizeof(uint32_t));
    uint32_t* tmp_off_arr = (uint32_t*) malloc(num_steps * sizeof(uint32_t));

    tbb::parallel_for( tbb::blocked_range<uint32_t>(0, num_steps, GRAIN_SIZE),
            [&](tbb::blocked_range<uint32_t> r){

            for (uint32_t i=r.begin(); i<r.end(); ++i){
                uint32_t index = GetKmerIndexAtPos(ref_str, start_addr+start_offset+(i*step), shape_size); 
                tmp_index_arr[i] = index;

                // valid index
                if (index != (uint32_t) INVALID_KMER) {
                    tmp_off_arr[i] = __sync_fetch_and_add( &index_table[index+1], 1);
                }
            }
        });

    InclusivePrefixScan(index_table, index_table_size);

    uint32_t num_index = index_table[index_table_size-1];
    
    pos_table = (uint32_t*) malloc(num_index * sizeof(uint32_t));

    tbb::parallel_for( tbb::blocked_range<uint32_t>(0, num_steps, GRAIN_SIZE),
            [&](tbb::blocked_range<uint32_t> r){

            for (uint32_t i=r.begin(); i<r.end(); ++i){
                uint32_t index = tmp_index_arr[i];

                // valid index
                if (index != (uint32_t) INVALID_KMER) {
                    uint32_t curr_idx = index_table[index] + tmp_off_arr[i]; 
                    pos_table[curr_idx] = start_offset+(i*step);
                }       
            }
        });

    SendSeedPosTable(index_table+1, index_table_size-1, pos_table, num_index);

    free(tmp_index_arr);
    free(tmp_off_arr);
    free(index_table);
    free(pos_table);
}
