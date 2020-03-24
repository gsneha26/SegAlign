#include "tbb/parallel_sort.h"
#include "seed_pos_table.h"

SeedPosTable::SeedPosTable() {
    ref_size_ = 0;
    kmer_size_ = 0;
    shape_size_ = 0;
}

int SeedPosTable::GetKmerSize() {
    return kmer_size_;
}

int SeedPosTable::GetShapeSize() {
    return shape_size_;
}

SeedPosTable::SeedPosTable(char* ref_str, size_t start_addr, uint32_t ref_length, std::string shape, uint32_t step) {

    int max_pos = 0;
    shape_size_ = shape.length(); 
    int kmer_size = 0;
    for (int i = 0; i < shape_size_; i++) {
        kmer_size += ((shape[i] == '1') || (shape[i] == 'T'));
    }
    
    assert(kmer_size <= 15);
    assert(kmer_size > 3); 

    kmer_size_ = kmer_size;
    ref_size_ = ref_length;

    GenerateShapePos(shape);

    uint32_t pos_table_size = ref_size_ - kmer_size_;
    assert(pos_table_size < ((uint64_t)1 << 32));

    index_table_size_ = ((uint32_t)1 << 2*kmer_size) + 1;
    index_table_ = (uint32_t*) calloc(index_table_size_, sizeof(uint32_t));
    pos_table_ = (uint64_t*) calloc(pos_table_size, sizeof(uint64_t));

    uint32_t num_index = 0;

    for (uint32_t i = 0; i < pos_table_size; i+=step) {
        uint32_t index = GetKmerIndexAtPos(ref_str, start_addr+i, shape_size_); 

        // valid index
        if (index != (1 << 31)) {
            pos_table_[num_index++] = ((uint64_t)index << 32) + i;
        }
    }

    tbb::parallel_sort(pos_table_, pos_table_+num_index);

    tmp_pos_table_ = (uint32_t*) calloc(num_index, sizeof(uint32_t));
    uint32_t curr_index = 0;
    uint32_t seed, pos; 
    uint32_t prev_i = 0;
    uint32_t num_i = 0;

    for (uint32_t i = 0; i < num_index; i++) {
        pos  = ((pos_table_[i] << 32) >> 32);
        seed = (pos_table_[i] >> 32);
        pos_table_[i] = pos;
        tmp_pos_table_[i] = pos;
        if (seed > curr_index) {
            for (uint32_t s = curr_index; s < seed; s++) {
                index_table_[s] = i;
                num_i = i - prev_i;
                if(num_i > max_pos)
                    max_pos = num_i;
                
                prev_i = i;
            }
            curr_index = seed;
        }
    }
    
    for (uint32_t i = curr_index; i < index_table_size_; i++) {
        index_table_[i] = num_index;
    }

    g_SendSeedPosTable(index_table_, index_table_size_, tmp_pos_table_, num_index, max_pos);

    free(index_table_);
    free(pos_table_);
    free(tmp_pos_table_);
}

SeedPosTable::~SeedPosTable() {
    free(index_table_);
    free(pos_table_);
    free(tmp_pos_table_);
}
