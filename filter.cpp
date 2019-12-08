#include <atomic>
#include <iostream>
#include "graph.h"

std::atomic<uint64_t> filter_body::num_filter_tiles(0);
std::atomic<uint64_t> filter_body::num_anchors(0);

extender_input filter_body::operator()(filter_input input){

    auto &payload = get<0>(input);

    auto &read = get<0>(payload);

    auto &data = get<1>(payload);

    size_t token = get<1>(input);
    filter_output fwOutput;
    filter_output rcOutput;
    


    auto &fwHits = data.fwHits;
    
    const size_t max_requests = (1 << 18); 

    uint8_t align_fields = 0;

    for (size_t b = 0; b < fwHits.size(); b += max_requests) {

        size_t first_hit = b;
        size_t last_hit = std::min(b + max_requests, fwHits.size());

        size_t num_requests = last_hit - first_hit;
        if (num_requests > 0) {
            std::vector<filter_tile> tiles;
            tiles.clear();
            for (size_t c = first_hit, r = 0; c < last_hit; c++, r++){
                uint32_t hit    = fwHits[c].reference_offset;
                uint32_t offset = fwHits[c].query_offset;
                
                size_t chr_id = std::upper_bound(r_chr_coord.cbegin(), r_chr_coord.cend(), hit) - r_chr_coord.cbegin() - 1;
                uint32_t chr_start = r_chr_coord[chr_id];
                uint32_t chr_end = chr_start + r_chr_len[chr_id];
                
                const size_t read_len = read.seq.size();
                char *read_char = (char *)read.seq.data();

                uint32_t ref_tile_start = 0;
                uint32_t query_tile_start = 0;
                
                uint32_t ref_tile_size = 300;
                uint32_t query_tile_size = 300;

                size_t ref_offset = ref_tile_start;
                size_t query_offset = (read_char - g_DRAM->buffer) - g_DRAM->referenceSize + query_tile_start;
                
                filter_tile tile = filter_tile(ref_offset, query_offset, ref_tile_size, query_tile_size, query_tile_start);
                tiles.push_back(tile);
                num_filter_tiles++;

            }

            std::vector<tile_output> f_op = g_SendBatchRequest(tiles, align_fields, cfg.xdrop_threshold);
            for (auto a: f_op) {
                int batch_id = a.batch_id;
                int score = a.tile_score;
                uint32_t ro = tiles[batch_id].ref_offset + a.max_ref_offset;
                uint32_t qo = tiles[batch_id].query_tile_start + a.max_query_offset;
                fwOutput.push_back(anchor(ro, qo, score));
            }
            num_anchors += f_op.size();
        }
    }

    auto &rcHits = data.rcHits;
    align_fields = reverse_query + complement_query;
    
    for (size_t b = 0; b < rcHits.size(); b += max_requests) {
        size_t first_hit = b;
        size_t last_hit = std::min(b + max_requests, rcHits.size());

        size_t num_requests = last_hit - first_hit;
        if (num_requests > 0) {
            std::vector<filter_tile> tiles;
            tiles.clear();
            for (size_t c = first_hit, r = 0; c < last_hit; c++, r++){
                uint32_t hit    = rcHits[c].reference_offset;
                uint32_t offset = rcHits[c].query_offset;
                
                size_t chr_id = std::upper_bound(r_chr_coord.cbegin(), r_chr_coord.cend(), hit) - r_chr_coord.cbegin() - 1;
                uint32_t chr_start = r_chr_coord[chr_id];
                uint32_t chr_end = chr_start + r_chr_len[chr_id];
                
                const size_t read_len = read.seq.size();
                char *read_char = (char *)read.seq.data();

                uint32_t ref_tile_start = 0;
                uint32_t query_tile_start =0;
                
                uint32_t ref_tile_size = 300;
                uint32_t query_tile_size = 300;

                size_t ref_offset = ref_tile_start;
                size_t query_offset = (read_char - g_DRAM->buffer) - g_DRAM->referenceSize + (read_len - (query_tile_start + query_tile_size));

                filter_tile tile = filter_tile(ref_offset, query_offset, ref_tile_size, query_tile_size, query_tile_start);
                tiles.push_back(tile);
                num_filter_tiles++;

            }

            std::vector<tile_output> f_op = g_SendBatchRequest(tiles, align_fields, cfg.xdrop_threshold);
            for (auto a: f_op) {
                int batch_id = a.batch_id;
                int score = a.tile_score;
                uint32_t ro = tiles[batch_id].ref_offset + a.max_ref_offset;
                uint32_t qo = tiles[batch_id].query_tile_start + a.max_query_offset;
                rcOutput.push_back(anchor(ro, qo, score));
            }
            num_anchors += f_op.size();
        }
    }

    std::sort(fwOutput.begin(), fwOutput.end(), CompareAnchors);
    std::sort(rcOutput.begin(), rcOutput.end(), CompareAnchors);

    return extender_input(extender_payload(read, fwOutput, rcOutput), token);
}
