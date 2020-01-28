#include "graph.h"
#include "ntcoding.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <atomic>

#include "tbb/parallel_for_each.h"

std::atomic<uint64_t> seeder_body::num_seed_hits(0);
std::atomic<uint64_t> seeder_body::num_seeds(0);
std::atomic<uint64_t> seeder_body::num_hsps(0);

printer_input seeder_body::operator()(seeder_input input) {

    seeder_payload &payload = get<0>(input);

    auto &query_chrom = get<0>(payload);
    
    auto &data = get<1>(payload);

    size_t token = get<1>(input);

    std::vector<hsp> fw_segments;
    std::vector<hsp> rc_segments;
    fw_segments.clear();
    rc_segments.clear();

    uint64_t index = 0;
    uint64_t transition_index = 0;
    uint64_t seed_offset;
    uint32_t start_pos = data.start;
    uint32_t end_pos = data.end;
    uint32_t num_invoked = data.num_invoked;
    uint32_t num_intervals = data.num_intervals;

    char* query = (char*) query_chrom.seq.data();
                    
    fprintf (stderr, "Chromosome %s interval %u/%u (%u:%u) \n", query_chrom.description.c_str(), num_invoked, num_intervals, start_pos, end_pos);

    for (uint32_t i = start_pos; i < end_pos; i += cfg.chunk_size) {

        //end position
        uint32_t e = std::min(i + cfg.chunk_size, end_pos);
        std::vector<uint64_t> seed_offset_vector;
        seed_offset_vector.clear();

        //start to end position in the chunk
        for (uint32_t j = i; j < e; j++) {
            index = GetKmerIndexAtPos(query, j);
            if (index != ((uint32_t) 1 << 31)) {
                seed_offset = (index << 32) + j;
                seed_offset_vector.push_back(seed_offset); 

                if (cfg.use_transition) {
                    for (int t=0; t < sa->GetKmerSize(); t++) {
                        if (IsTransitionAtPos(t) == 1) {
                            transition_index = (index ^ (TRANSITION_MASK << (2*t)));
                            seed_offset = (transition_index << 32) + j;
                            seed_offset_vector.push_back(seed_offset); 
                        }
                    }
                }
            }
        }
        
        fw_segments = g_SeedAndFilter(seed_offset_vector, false);
    }

    char* rc_query = (char*) query_chrom.rc_seq.data();
    for (uint32_t i = start_pos; i < end_pos; i += cfg.chunk_size) {
        uint32_t e = std::min(i + cfg.chunk_size, end_pos);
        std::vector<uint64_t> seed_offset_vector;
        seed_offset_vector.clear();
        for (uint32_t j = i; j < e; j++) {
            index = GetKmerIndexAtPos(rc_query, j);
            if (index != ((uint32_t) 1 << 31)) {
                seed_offset = (index << 32) + j;
                seed_offset_vector.push_back(seed_offset); 
                if (cfg.use_transition) {
                    for (int t=0; t < sa->GetKmerSize(); t++) {
                        if (IsTransitionAtPos(t) == 1) {
                            transition_index = (index ^ (TRANSITION_MASK << (2*t)));
                            seed_offset = (transition_index << 32) + j;
                            seed_offset_vector.push_back(seed_offset); 
                        }
                    }
                }
            }
        }

        rc_segments = g_SeedAndFilter(seed_offset_vector, true);
    }

	return printer_input(printer_payload(fw_segments, rc_segments), token);
}

