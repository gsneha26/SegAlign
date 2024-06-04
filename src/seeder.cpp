#include "graph.h"
#include "ntcoding.h"
#include "seed_filter.h"
#include "store.h"

std::atomic<uint64_t> seeder_body::num_seed_hits(0);
std::atomic<uint64_t> seeder_body::num_seeds(0);
std::atomic<uint64_t> seeder_body::num_hsps(0);
std::atomic<uint32_t> seeder_body::total_xdrop(0);
std::atomic<uint32_t> seeder_body::num_seeded_regions[BUFFER_DEPTH]={};

printer_input seeder_body::operator()(seeder_input input) {

    auto &payload = std::get<0>(input);
    size_t token  = std::get<1>(input);

    auto &block_data    = std::get<0>(payload);
    auto &interval_data = std::get<1>(payload);

    int r_block_index    = block_data.r_index;
    int q_block_index    = block_data.q_index;
    size_t r_block_start = block_data.r_start;
    size_t q_block_start = block_data.q_start;
    uint32_t r_block_len = block_data.r_len;
    uint32_t q_block_len = block_data.q_len;

    uint32_t q_inter_start = interval_data.start;
    uint32_t q_inter_end   = interval_data.end;
    uint32_t buffer        = interval_data.buffer;
    uint32_t num_invoked   = interval_data.num_invoked;
    uint32_t num_intervals = interval_data.num_intervals;

    uint32_t rc_q_inter_start = q_block_len - q_inter_end;
    uint32_t rc_q_inter_end   = q_block_len - q_inter_start;

    uint64_t kmer_index;
    uint64_t transition_index;
    uint64_t seed_offset;

    std::vector<segmentPair> fw_hsps;
    std::vector<segmentPair> rc_hsps;
    fw_hsps.clear();
    rc_hsps.clear();

    fprintf (stderr, "Query block %u, interval %u/%u (%u:%u) with buffer %u\n", q_block_index, num_invoked, num_intervals, q_inter_start, q_inter_end, buffer);

    if(cfg.strand == "plus" || cfg.strand == "both"){
        for (uint32_t i = q_inter_start; i < q_inter_end; i += cfg.wga_chunk_size) {

            //end position
            uint32_t e = std::min(i + cfg.wga_chunk_size, q_inter_end);

            std::vector<uint64_t> seed_offset_vector;
            seed_offset_vector.clear();

            //start to end position in the chunk
            for (uint32_t j = i; j < e; j++) {

                kmer_index = GetKmerIndexAtPos(query_DRAM->buffer, q_block_start+j, cfg.seed.size);
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
                std::vector<segmentPair> anchors = g_SeedAndFilter(seed_offset_vector, false, buffer);
                seeder_body::num_seed_hits += anchors[0].score;
                if(anchors.size() > 1){
                    fw_hsps.insert(fw_hsps.end(), anchors.begin()+1, anchors.end());
                    seeder_body::num_hsps += anchors.size()-1;
                }
            }
        }
    }

    if(cfg.strand == "minus" || cfg.strand == "both"){
        for (uint32_t i = rc_q_inter_start; i < rc_q_inter_end; i += cfg.wga_chunk_size) {
            uint32_t e = std::min(i + cfg.wga_chunk_size, rc_q_inter_end);

            std::vector<uint64_t> seed_offset_vector;
            seed_offset_vector.clear();
            for (uint32_t j = i; j < e; j++) {
                kmer_index = GetKmerIndexAtPos(query_rc_DRAM->buffer, q_block_start+j, cfg.seed.size);
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
                std::vector<segmentPair> anchors = g_SeedAndFilter(seed_offset_vector, true, buffer);
                seeder_body::num_seed_hits += anchors[0].score;
                if(anchors.size() > 1){
                    rc_hsps.insert(rc_hsps.end(), anchors.begin()+1, anchors.end());
                    seeder_body::num_hsps += anchors.size()-1;
                }
            }
        }
    }

    seeder_body::num_seeded_regions[buffer] += 1;
    seeder_body::total_xdrop += 1;

    return printer_input(printer_payload(payload, fw_hsps, rc_hsps), token);
}
