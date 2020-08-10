#include "graph.h"
#include "store.h"

std::atomic<uint64_t> seeder_body::num_seed_hits(0);
std::atomic<uint64_t> seeder_body::num_seeds(0);
std::atomic<uint64_t> seeder_body::num_hsps(0);
std::atomic<uint32_t> seeder_body::total_xdrop(0);
std::atomic<uint32_t> seeder_body::num_seeded_regions(0);

printer_input seeder_body::operator()(seeder_input input) {

    seeder_payload &payload = get<0>(input);

    auto &chrom = get<0>(payload);
    
    auto &data = get<1>(payload);

    size_t token = get<1>(input);

    std::vector<segment> fw_segments;
    std::vector<segment> rc_segments;
    fw_segments.clear();
    rc_segments.clear();

    uint64_t index = 0;
    uint64_t transition_index = 0;
    uint64_t seed_offset;
    size_t block_start = chrom.r_start;
    uint32_t block_len = chrom.r_len;
    uint32_t block_index = chrom.r_block_index;
    uint32_t start_pos = data.start;
    uint32_t end_pos = data.end;
    uint32_t ref_start = data.ref_start;
    uint32_t ref_end = data.ref_end;
    uint32_t num_invoked = data.num_invoked;
    uint32_t num_intervals = data.num_intervals;
    uint32_t start_pos_rc = block_len - 1 - end_pos;
    uint32_t end_pos_rc = block_len - 1 - start_pos;
    uint32_t rc_block_start = cfg.ref_len - 1 - block_start - (block_len - 1);

    fprintf (stderr, "Chromosome block %u interval %u/%u (%u:%u) with ref (%u:%u) rc (%u:%u)\n", block_index, num_invoked, num_intervals, block_start+start_pos, block_start+end_pos, ref_start, ref_end, rc_block_start+start_pos_rc, rc_block_start+end_pos_rc);
//    fprintf (stderr, "Chromosome block %u interval %u/%u (%u:%u) with ref (%u:%u) rc (%u:%u)\n%u %u %u %u %u %u %u\n", block_index, num_invoked, num_intervals, block_start+start_pos, block_start+end_pos, ref_start, ref_end, rc_block_start+start_pos_rc, rc_block_start+end_pos_rc, block_start, block_len, start_pos, end_pos, start_pos_rc, end_pos_rc, rc_block_start);

    if(cfg.strand == "plus" || cfg.strand == "both"){
        for (uint32_t i = start_pos; i < end_pos; i += cfg.wga_chunk_size) {

            //end position
            uint32_t e = std::min(i + cfg.wga_chunk_size, end_pos);

            std::vector<uint64_t> seed_offset_vector;
            seed_offset_vector.clear();

            //start to end position in the chunk
            for (uint32_t j = i; j < e; j++) {

                index = GetKmerIndexAtPos(seq_DRAM->buffer, block_start+j, cfg.seed_size);
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
                std::vector<segment> anchors = g_SeedAndFilter(seed_offset_vector, false, ref_start, ref_end);
                seeder_body::num_seed_hits += anchors[0].score;
                if(anchors.size() > 1){
                    fw_segments.insert(fw_segments.end(), anchors.begin()+1, anchors.end());
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
                index = GetKmerIndexAtPos(seq_rc_DRAM->buffer, (rc_block_start+j), cfg.seed_size);
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
                std::vector<segment> anchors = g_SeedAndFilter(seed_offset_vector, true, ref_start, ref_end);
                seeder_body::num_seed_hits += anchors[0].score;
                if(anchors.size() > 1){
                    rc_segments.insert(rc_segments.end(), anchors.begin()+1, anchors.end());
                    seeder_body::num_hsps += anchors.size()-1;
                }
            }
        }
    }

    seeder_body::num_seeded_regions += 1;
    seeder_body::total_xdrop += 1;

    return printer_input(printer_payload(chrom, num_invoked, fw_segments, rc_segments), token);
}
