#include <algorithm>
#include <iterator>
#include <string>
#include <map>
#include <mutex>
#include "graph.h"
#include "store.h"

std::mutex io_lock;

void segment_printer_body::operator()(printer_input input, printer_node::output_ports_type & op){

    auto &payload = get<0>(input); 
    size_t token  = get<1>(input);
 
    auto &block_data  = get<0>(payload);
    auto &index       = get<1>(payload);
    auto &out_hsps = get<2>(payload);

    int block_index    = block_data.index;
    size_t block_start = block_data.start;
    uint32_t block_len = block_data.len;
    size_t block_end   = block_start + block_len;

    std::string segment_filename;
    std::string cmd;
    
    uint32_t num_hsps  = out_hsps.size();

    size_t q_index;
    size_t seg_q_start;

    if(num_hsps > 0){

        size_t start_chr = std::upper_bound(chr_start.cbegin(), chr_start.cend(), block_start) - chr_start.cbegin() - 1; 
        size_t end_chr = std::upper_bound(chr_start.cbegin(), chr_start.cend(), block_end) - chr_start.cbegin() - 1; 

        uint32_t curr_q_chr_index = start_chr;
        std::string curr_q_chr    = chr_name[start_chr];
        size_t curr_q_chr_start   = chr_start[start_chr]; 
        size_t curr_q_chr_end     = curr_q_chr_start + chr_len[start_chr];

        std::string out_str;
        FILE* segmentFile;

        segment_filename = "tmp"+std::to_string(index)+".block"+std::to_string(block_index)+".intervals"; 
        segmentFile = fopen(segment_filename.c_str(), "w");

        for (auto e: out_hsps) {
            seg_q_start = block_start + e.query_start;

            if(seg_q_start < curr_q_chr_start || seg_q_start >= curr_q_chr_end){
                q_index = std::upper_bound(chr_start.cbegin(), chr_start.cend(), seg_q_start) - chr_start.cbegin() - 1;
                curr_q_chr_index = q_index;
                curr_q_chr       = chr_name[curr_q_chr_index];
                curr_q_chr_start = chr_start[curr_q_chr_index];
                curr_q_chr_end   = curr_q_chr_start + chr_len[curr_q_chr_index];
            }

            fprintf(segmentFile, "%s\t%lu\t%lu\n", curr_q_chr.c_str(), seg_q_start-curr_q_chr_start, seg_q_start+e.len+1-curr_q_chr_start);
        }

        fclose(segmentFile);
    }

    get<0>(op).try_put(token);
};
