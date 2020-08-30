#include "graph.h"
#include "store.h"
#include "ntcoding.h"
#include "seed_filter.h"

std::atomic<uint64_t> seeder_body::num_seed_hits(0);
std::atomic<uint64_t> seeder_body::num_seeds(0);
std::atomic<uint64_t> seeder_body::num_hsps(0);
std::atomic<uint32_t> seeder_body::total_xdrop(0);
std::atomic<uint32_t> seeder_body::num_seeded_regions(0);

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

    std::vector<segmentPair> fw_hsps;
    std::vector<segmentPair> rc_hsps;
    fw_hsps.clear();
    rc_hsps.clear();

    fprintf (stderr, "Chromosome block %u interval %u/%u (%lu:%lu) with ref (%u:%u) rc (%lu:%lu)\n", block_index, num_invoked, num_intervals, block_start+start_pos, block_start+end_pos, ref_start, ref_end, rc_block_start+start_pos_rc, rc_block_start+end_pos_rc);

    if(cfg.strand == "plus" || cfg.strand == "both"){
        for (uint32_t i = start_pos; i < end_pos; i += cfg.wga_chunk_size) {

            //end position
            uint32_t e = std::min(i + cfg.wga_chunk_size, end_pos);

            std::vector<uint64_t> seed_offset_vector;
            seed_offset_vector.clear();

            //start to end position in the chunk
            for (uint32_t j = i; j < e; j++) {

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
                seeder_body::num_seed_hits += anchors[0].score;
                if(anchors.size() > 1){
                    fw_hsps.insert(fw_hsps.end(), anchors.begin()+1, anchors.end());
                    seeder_body::num_hsps += anchors.size()-1;
                }
            }
        }
    }

    if(cfg.strand == "minus" || cfg.strand == "both"){
        for (uint32_t i = start_pos_rc; i < end_pos_rc; i += cfg.wga_chunk_size) {
            uint32_t e = std::min(i + cfg.wga_chunk_size, end_pos_rc);

            std::vector<uint64_t> seed_offset_vector;
            seed_offset_vector.clear();
            for (uint32_t j = i; j < e; j++) {
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
                seeder_body::num_seed_hits += anchors[0].score;
                if(anchors.size() > 1){
                    rc_hsps.insert(rc_hsps.end(), anchors.begin()+1, anchors.end());
                    seeder_body::num_hsps += anchors.size()-1;
                }
            }
        }
    }

    seeder_body::num_seeded_regions += 1;
    seeder_body::total_xdrop += 1;

    return printer_input(printer_payload(block_data, num_invoked, fw_hsps, rc_hsps), token);
}
