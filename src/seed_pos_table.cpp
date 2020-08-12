#include "tbb/parallel_sort.h"
#include <tbb/mutex.h>
#include "seed_pos_table.h"

tbb::mutex l1;

SeedPosTable::SeedPosTable(){
    ref_size_ = 0;
    kmer_size_ = 0;
    shape_size_ = 0;
}

int SeedPosTable::GetKmerSize(){
    return kmer_size_;
}

SeedPosTable::SeedPosTable(char* ref_str, size_t start_addr, uint32_t ref_length, std::string shape, uint32_t step) {

    shape_size_ = shape.length(); 
    int kmer_size = 0;
    for (int i = 0; i < shape_size_; i++){
        kmer_size += ((shape[i] == '1') || (shape[i] == 'T'));
    }
    
    assert(kmer_size <= 15);
    assert(kmer_size > 3); 

    kmer_size_ = kmer_size;
    ref_size_ = ref_length;

    GenerateShapePos(shape);

    uint32_t offset = (shape_size_+1)%step;
    uint32_t start_offset = step - offset;
    
    index_table_size_ = ((uint32_t)1 << 2*kmer_size) + 1;
    index_table_ = (uint32_t*) calloc(index_table_size_, sizeof(uint32_t));

    uint32_t num_steps = (ref_length - shape_size_ + offset) / step;

    uint32_t* tmp_index_arr = (uint32_t*) malloc(num_steps * sizeof(uint32_t));
    uint32_t* tmp_off_arr = (uint32_t*) malloc(num_steps * sizeof(uint32_t));

    tbb::parallel_for( tbb::blocked_range<uint32_t>(0, num_steps, GRAIN_SIZE),
            [&](tbb::blocked_range<uint32_t> r){

            for (uint32_t i=r.begin(); i<r.end(); ++i){
                uint32_t index = GetKmerIndexAtPos(ref_str, start_addr+start_offset+(i*step), shape_size_); 
                tmp_index_arr[i] = index;

                // valid index
                if (index != INVALID_KMER) {
                    tmp_off_arr[i] = __sync_fetch_and_add( &index_table_[index+1], 1);
                }
            }
        });

    g_InclusivePrefixScan(index_table_, index_table_size_);

    uint32_t num_index = index_table_[index_table_size_-1];
    
    pos_table_ = (uint32_t*) malloc(num_index * sizeof(uint32_t));

    tbb::parallel_for( tbb::blocked_range<uint32_t>(0, num_steps, GRAIN_SIZE),
            [&](tbb::blocked_range<uint32_t> r){

            for (uint32_t i=r.begin(); i<r.end(); ++i){
                uint32_t index = tmp_index_arr[i];

                // valid index
                if (index != INVALID_KMER) {
                    uint32_t curr_idx = index_table_[index] + tmp_off_arr[i]; 
                    pos_table_[curr_idx] = start_offset+(i*step);
                }       
            }
        });

    g_SendSeedPosTable(index_table_+1, index_table_size_-1, pos_table_, num_index);

    free(tmp_index_arr);
    free(tmp_off_arr);
    free(index_table_);
    free(pos_table_);
}

SeedPosTable::~SeedPosTable() {
    free(index_table_);
    free(pos_table_);
}
