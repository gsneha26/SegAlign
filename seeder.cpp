#include "graph.h"
#include "store.h"
#include <atomic>

std::atomic<uint64_t> seeder_body::num_seed_hits(0);
std::atomic<uint64_t> seeder_body::num_seeds(0);
std::atomic<uint64_t> seeder_body::num_hsps(0);
std::atomic<uint32_t> seeder_body::total_xdrop(0);
std::atomic<uint32_t> seeder_body::num_seeded_regions[BUFFER_DEPTH]={};

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
    uint32_t q_len = chrom.q_len;
    size_t q_start = chrom.q_start;
    size_t r_start = chrom.r_start;
    uint32_t q_block_index = chrom.block_index;

//    fprintf (stderr, "Chromosome block %u interval %u/%u (%u:%u) with buffer %u\n", q_block_index, num_invoked, num_intervals, start_pos, end_pos, buffer);

    for (uint32_t i = start_pos; i < end_pos; i += cfg.wga_chunk_size) {

        //end position
        uint32_t e = std::min(i + cfg.wga_chunk_size, end_pos);

        std::vector<uint64_t> seed_offset_vector;
        seed_offset_vector.clear();

        //start to end position in the chunk
        for (uint32_t j = i; j < e; j++) {
            index = GetKmerIndexAtPos(query_DRAM->buffer, q_start+j, cfg.seed_size);
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
            seeder_body::num_seeds += seed_offset_vector.size();
            std::vector<hsp> anchors = g_SeedAndFilter(seed_offset_vector, false, buffer, cfg.seed_size, cfg.xdrop, cfg.hspthresh, cfg.noentropy, cfg.nounique); 
            seeder_body::num_seed_hits += anchors[0].score;
            if(anchors.size() > 1){
                fw_segments.insert(fw_segments.end(), anchors.begin()+1, anchors.end());
                seeder_body::num_hsps += anchors.size()-1;
            }
        }
    }

    for (uint32_t i = q_len - end_pos; i < q_len - start_pos; i += cfg.wga_chunk_size) {
        uint32_t e = std::min(i + cfg.wga_chunk_size, q_len - start_pos);

        std::vector<uint64_t> seed_offset_vector;
        seed_offset_vector.clear();
        for (uint32_t j = i; j < e; j++) {
            index = GetKmerIndexAtPos(query_rc_DRAM->buffer, q_start+j, cfg.seed_size);
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
            seeder_body::num_seeds += seed_offset_vector.size();
            std::vector<hsp> anchors = g_SeedAndFilter(seed_offset_vector, true, buffer, cfg.seed_size, cfg.xdrop, cfg.hspthresh, cfg.noentropy, cfg.nounique); 
            seeder_body::num_seed_hits += anchors[0].score;
            if(anchors.size() > 1){
                rc_segments.insert(rc_segments.end(), anchors.begin()+1, anchors.end());
                seeder_body::num_hsps += anchors.size()-1;
            }
        }
    }

    seeder_body::num_seeded_regions[buffer] += 1;
    seeder_body::total_xdrop += 1;

    return printer_input(printer_payload(num_invoked, fw_segments, rc_segments, q_block_index, r_start, q_start), token);
}
