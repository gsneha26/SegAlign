#include "ntcoding.h"
#include "seed_pos_table.h"
#include "GPU.h"
#include "tbb/parallel_scan.h"
#include "tbb/parallel_sort.h"
#include "tbb/blocked_range.h"
#include "tbb/scalable_allocator.h"
#include <algorithm>
#include <atomic>

SeedPosTable::SeedPosTable() {
    ref_size_ = 0;
    kmer_size_ = 0;
    shape_size_ = 0;
    bin_size_ = 0;
}

int SeedPosTable::GetKmerSize() {
    return kmer_size_;
}

int SeedPosTable::GetShapeSize() {
    return shape_size_;
}

SeedPosTable::SeedPosTable(char* ref_str, uint32_t ref_length, std::string shape, int bin_size) {
    shape_size_ = shape.length(); 
    int kmer_size = 0;
    for (int i = 0; i < shape_size_; i++) {
        kmer_size += ((shape[i] == '1') || (shape[i] == 'T'));
    }
    
    assert(kmer_size <= 15);
    assert(kmer_size > 3); 

    kmer_size_ = kmer_size;
    ref_size_ = ref_length;
    bin_size_ = bin_size;

    GenerateShapePos(shape);

    uint32_t pos_table_size = ref_size_ - kmer_size_;
    assert(pos_table_size < ((uint64_t)1 << 32));

    index_table_size_ = ((uint32_t)1 << 2*kmer_size) + 1;
    index_table_ = (uint32_t*) calloc(index_table_size_, sizeof(uint32_t));
    pos_table_ = (uint64_t*) calloc(pos_table_size, sizeof(uint64_t));

    uint32_t num_index = 0;

    for (uint32_t i = 0; i < pos_table_size; i++) {
        uint32_t index = GetKmerIndexAtPos(ref_str, i); 

        // valid index
        if (index != (1 << 31)) {
            pos_table_[num_index++] = ((uint64_t)index << 32) + i;
        }
    }


    tbb::parallel_sort(pos_table_, pos_table_+num_index);

    uint32_t curr_index = 0;
    uint32_t seed, pos; 

    for (uint32_t i = 0; i < num_index; i++) {
        pos  = ((pos_table_[i] << 32) >> 32);
        seed = (pos_table_[i] >> 32);
        pos_table_[i] = pos;
        if (seed > curr_index) {
            for (uint32_t s = curr_index; s < seed; s++) {
                index_table_[s] = i;
            }
            curr_index = seed;
        }
    }
    
    for (uint32_t i = curr_index; i < index_table_size_; i++) {
        index_table_[i] = num_index;
    }

    g_SendSeedPosTable(index_table_, index_table_size_, pos_table_, num_index);
}

SeedPosTable::~SeedPosTable() {
    free(index_table_);
    free(pos_table_);
}


std::vector<seed_hit> SeedPosTable::DSOFT(std::vector<uint64_t> seed_offset_vector, int threshold, uint32_t chunk_offset) {
    std::vector<seed_hit> seed_hits;
    seed_hits.clear();
    
//	std::vector<Hits> hits_array;
//    hits_array.clear();
//
    uint32_t num_seeds = seed_offset_vector.size();
//
//    for (uint32_t i = 0; i < num_seeds; i++) {
//		uint32_t index  = (seed_offset_vector[i] >> 32);
//		uint32_t offset = ((seed_offset_vector[i] << 32) >> 32);
//
//		uint32_t start_index = (index == 0) ?  0 : index_table_[index-1];
//		uint32_t end_index = index_table_[index];
//        
//        for (uint32_t j = start_index; j < end_index; j++) {
//            uint32_t hit = pos_table_[j];
//            seed_hit sh;
//            sh.reference_offset = hit;
//            sh.query_offset = offset;
//            seed_hits.push_back(sh);
////            if (hit >= offset) {
////                uint32_t bin = ((hit - offset) / bin_size_);
////                uint64_t bin_offset = ((uint64_t)bin << 32) + offset;
////                hits_array.push_back(Hits(bin_offset, hit));
////            }
//        }
//    }
//	
//    std::stable_sort(hits_array.begin(), hits_array.end(), CompareHits);
//	
//    uint32_t num_hits = hits_array.size();
//	
//    uint32_t last_bin = (1 << 31);
//	uint32_t last_offset = 0;
//	uint32_t curr_count = 0;
//
//	for (int i = 0; i < num_hits; i++) {
//		Hits next_hit = hits_array[i];
//		uint32_t offset = ((next_hit.bin_offset << 32) >> 32);
//		uint32_t bin = (next_hit.bin_offset >> 32);
//		uint32_t hit = next_hit.hit;
//		if (bin == last_bin) {
//            if (curr_count < threshold) {
//                curr_count = ((offset - last_offset > kmer_size_) || (curr_count == 0)) ? curr_count + kmer_size_ : curr_count + (offset - last_offset);
//                if (curr_count >= threshold) {
//                    seed_hit sh;
//                    sh.reference_offset = hit;
//                    sh.query_offset = offset + chunk_offset;
//                    seed_hits.push_back(sh);
//                }
//            }
//        }
//        else {
//            last_bin = bin;
//            curr_count = kmer_size_;
//            if (curr_count >= threshold) {
//                seed_hit sh;
//                sh.reference_offset = hit;
//                sh.query_offset = offset + chunk_offset;
//                seed_hits.push_back(sh);
//            }
//        }
//		last_offset = offset;
//	}


    return seed_hits;
}

