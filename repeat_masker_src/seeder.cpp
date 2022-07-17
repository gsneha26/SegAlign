#include "graph.h"
#include <algorithm>
#include <assert.h>
#include "store.h"
#include "ntcoding.h"
#include "seed_filter.h"
#include "tbb/scalable_allocator.h"

std::atomic<uint64_t> seeder_body::num_seed_hits(0);
std::atomic<uint64_t> seeder_body::num_seeds(0);
std::atomic<uint64_t> seeder_body::num_hsps(0);
std::atomic<uint32_t> seeder_body::total_xdrop(0);
std::atomic<uint32_t> seeder_body::num_seeded_regions(0);

bool sort_hsp(segmentPair x, segmentPair y){

    if(x.query_start < y.query_start)
        return true;
    else if(x.query_start == y.query_start){
        if(x.len < y.len)
            return true;
        else
            return false;
    }
    return false;
}

printer_input seeder_body::operator()(seeder_input input) {

    auto &payload = get<0>(input);
    size_t token  = get<1>(input);

    auto &block_data    = get<0>(payload);
    auto &interval_data = get<1>(payload);

    size_t block_start = block_data.start;
    uint32_t block_len = block_data.len;
    int block_index    = block_data.index;

    uint32_t start_pos = interval_data.start;
    uint32_t end_pos   = interval_data.end;
    uint32_t ref_start = interval_data.ref_start;
    uint32_t ref_end   = interval_data.ref_end;
    uint32_t num_invoked   = interval_data.num_invoked;
    uint32_t num_intervals = interval_data.num_intervals;

    uint32_t start_pos_rc = block_len - 1 - end_pos;
    uint32_t end_pos_rc   = block_len - 1 - start_pos;
    size_t rc_block_start = cfg.seq_len - 1 - block_start - (block_len - 1);

    uint64_t kmer_index;
    uint64_t transition_index;
    uint64_t seed_offset;

    std::vector<segmentPair> total_hsps;
    total_hsps.clear();

    uint8_t* int_count;
    int_count = (uint8_t*) scalable_calloc(block_len, sizeof(uint8_t)); 

    std::vector<Segment> total_intervals;
    total_intervals.clear();

    int32_t start;
    int32_t end;

    uint32_t old_num_hsps = 0;
    uint32_t new_num_hsps = 0;

    fprintf (stderr, "Chromosome block %u interval %u/%u (%lu:%lu) with ref (%u:%u) rc (%lu:%lu)\n", block_index, num_invoked, num_intervals, block_start+start_pos, block_start+end_pos, ref_start, ref_end, rc_block_start+start_pos_rc, rc_block_start+end_pos_rc);

    for (uint32_t i = start_pos; i < end_pos; i += cfg.wga_chunk_size) {

        //chunk limit positions
        start = i;
        end = std::min(start + cfg.wga_chunk_size, end_pos);

        if(cfg.strand == "plus" || cfg.strand == "both"){

            std::vector<uint64_t> seed_offset_vector;
            seed_offset_vector.clear();

            //start to end position in the chunk
            for (uint32_t j = start; j < end; j++) {

                kmer_index = GetKmerIndexAtPos(seq_DRAM->buffer, block_start+j, cfg.seed.size);
                if (kmer_index != ((uint32_t) 1 << 31)) {
                    seed_offset = (kmer_index << 32) + j;
                    seed_offset_vector.push_back(seed_offset); 

                    if (cfg.seed.transition) {
                        for (int t=0; t < cfg.seed.kmer_size; t++) {
                            if (IsTransitionAtPos(t) == 1) {
                                transition_index = (kmer_index ^ (TRANSITION_MASK << (2*t)));
                                seed_offset = (transition_index << 32) + j;
                                seed_offset_vector.push_back(seed_offset); 
                            }
                        }
                    }
                }
            }

            if(seed_offset_vector.size() > 0){
                seeder_body::num_seeds += seed_offset_vector.size();
                std::vector<segmentPair> anchors = g_SeedAndFilter(seed_offset_vector, false, ref_start, ref_end);
		uint64_t curr_hits = ((uint64_t)anchors[0].query_start << 32) + anchors[0].ref_start; 
                seeder_body::num_seed_hits += curr_hits;
                if(anchors.size() > 1){
                    total_hsps.insert(total_hsps.end(), anchors.begin()+1, anchors.end());
                    seeder_body::num_hsps += anchors.size()-1;
                }
            }
        }

        if(cfg.strand == "minus" || cfg.strand == "both"){

            //chunk limit positions
            start = block_len - 1 - end;
            end = std::min(start + cfg.wga_chunk_size, end_pos_rc);

            std::vector<uint64_t> seed_offset_vector;
            seed_offset_vector.clear();

            for (uint32_t j = start; j < end; j++) {
                kmer_index = GetKmerIndexAtPos(seq_rc_DRAM->buffer, (rc_block_start+j), cfg.seed.size);
                if (kmer_index != ((uint32_t) 1 << 31)) {
                    seed_offset = (kmer_index << 32) + j;
                    seed_offset_vector.push_back(seed_offset); 
                    if (cfg.seed.transition) {
                        for (int t=0; t < cfg.seed.kmer_size; t++) {
                            if (IsTransitionAtPos(t) == 1) {
                                transition_index = (kmer_index ^ (TRANSITION_MASK << (2*t)));
                                seed_offset = (transition_index << 32) + j;
                                seed_offset_vector.push_back(seed_offset); 
                            }
                        }
                    }
                }
            }

            if(seed_offset_vector.size() > 0){
                seeder_body::num_seeds += seed_offset_vector.size();
                std::vector<segmentPair> anchors = g_SeedAndFilter(seed_offset_vector, true, ref_start, ref_end);
		uint64_t curr_hits = ((uint64_t)anchors[0].query_start << 32) + anchors[0].ref_start; 
                seeder_body::num_seed_hits += curr_hits;
                if(anchors.size() > 1){
                    total_hsps.insert(total_hsps.end(), anchors.rbegin(), anchors.rend()-1);
                    seeder_body::num_hsps += anchors.size()-1;
                }
            }
        }

        if(total_hsps.size() > 0){

            std::sort(total_hsps.begin(), total_hsps.end(), sort_hsp);
	    for (auto hsp: total_hsps) {
                for(int j = hsp.query_start; j < (hsp.query_start + hsp.len); j++){
                    int_count[j]++;
                }
            }
            total_hsps.clear();
        }
    }

    int run = 0;
    int query_start = 0;
    int len = 0;

    for(int i = 0; i < block_len; i++){
        if(int_count[i] >= cfg.M){
            if(run == 0){
                run = 1;
                query_start = i;
            }
            len++;
        }
        else{
            if(run == 1){
                run = 0;
                Segment s;
                s.query_start = query_start;
                s.len = len;
                total_intervals.push_back(s);
            }
            query_start = 0;
            len = 0;
        }
    }

    scalable_free(int_count);

    seeder_body::num_seeded_regions += 1;
    seeder_body::total_xdrop += 1;

    return printer_input(printer_payload(block_data, num_invoked, total_intervals), token);
}
