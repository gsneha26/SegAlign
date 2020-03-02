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
std::atomic<uint32_t> seeder_body::num_seeded_regions0(0);
std::atomic<uint32_t> seeder_body::num_seeded_regions1(0);

printer_input seeder_body::operator()(seeder_input input) {

    seeder_payload &payload = get<0>(input);

    auto &chrom = get<0>(payload);
    
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
    uint32_t buffer = data.buffer;

    fprintf (stderr, "Chromosome %s interval %u/%u (%u:%u) with thread %d\n", chrom.query_chr.c_str(), num_invoked, num_intervals, start_pos, end_pos, token);

    char* query = (char*) chrom.q_seq.data();

    for (uint32_t i = start_pos; i < end_pos; i += WGA_CHUNK) {

        //end position
        uint32_t e = std::min(i + WGA_CHUNK, end_pos);
        std::vector<uint64_t> seed_offset_vector;
        seed_offset_vector.clear();

        //start to end position in the chunk
        for (uint32_t j = i; j < e; j++) {
            index = GetKmerIndexAtPos(query, j);
            if (index != ((uint32_t) 1 << 31)) {
                seed_offset = (index << 32) + j;
                seed_offset_vector.push_back(seed_offset); 

                if (cfg.transition) {
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

        if(seed_offset_vector.size() > 0){
            std::vector<hsp> anchors = g_SeedAndFilter(seed_offset_vector, false, buffer); 
            fw_segments.insert(fw_segments.end(), anchors.begin(), anchors.end());
            seeder_body::num_seeds += seed_offset_vector.size();
            seeder_body::num_hsps += anchors.size();
        }
    }

    char* rc_query = (char*) chrom.q_rc_seq.data();
    for (uint32_t i = start_pos; i < end_pos; i += WGA_CHUNK) {
        uint32_t e = std::min(i + WGA_CHUNK, end_pos);
        std::vector<uint64_t> seed_offset_vector;
        seed_offset_vector.clear();
        for (uint32_t j = i; j < e; j++) {
            index = GetKmerIndexAtPos(rc_query, j);
            if (index != ((uint32_t) 1 << 31)) {
                seed_offset = (index << 32) + j;
                seed_offset_vector.push_back(seed_offset); 
                if (cfg.transition) {
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

        if(seed_offset_vector.size() > 0){
            std::vector<hsp> anchors = g_SeedAndFilter(seed_offset_vector, true, buffer); 
            rc_segments.insert(rc_segments.end(), anchors.begin(), anchors.end());
            seeder_body::num_seeds += seed_offset_vector.size();
            seeder_body::num_hsps += anchors.size();
        }
    }

    if(buffer == 0){
        seeder_body::num_seeded_regions0 += 1;
    }
    else{
        seeder_body::num_seeded_regions1 += 1;
    }

    return printer_input(printer_payload(num_invoked, fw_segments, rc_segments, chrom.query_chr, chrom.ref_chr), token);
}

