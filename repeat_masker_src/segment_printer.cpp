#include <algorithm>
#include <iterator>
#include <string>
#include <map>
#include "graph.h"
#include "store.h"

void interval_printer_body::operator()(printer_input input, printer_node::output_ports_type & op){

    auto &payload = std::get<0>(input); 
    size_t token  = std::get<1>(input);
 
    auto &block_data  = std::get<0>(payload);
    auto &index       = std::get<1>(payload);
    auto &out_intervals    = std::get<2>(payload);

    int block_index    = block_data.index;
    size_t block_start = block_data.start;
    uint32_t block_len = block_data.len;
    size_t block_end   = block_start + block_len;

    std::string interval_filename;
    
    uint32_t num_intervals  = out_intervals.size();

    size_t q_index;
    size_t seg_q_start;

    if(num_intervals > 0){

        size_t start_chr = std::upper_bound(chr_start.cbegin(), chr_start.cend(), block_start) - chr_start.cbegin() - 1; 
        size_t end_chr = std::upper_bound(chr_start.cbegin(), chr_start.cend(), block_end) - chr_start.cbegin() - 1; 

        uint32_t curr_q_chr_index = start_chr;
        std::string curr_q_chr    = chr_name[start_chr];
        size_t curr_q_chr_start   = chr_start[start_chr]; 
        size_t curr_q_chr_end     = curr_q_chr_start + chr_len[start_chr];

        FILE* intervalFile;

        interval_filename = "tmp"+std::to_string(index)+".block"+std::to_string(block_index)+".intervals"; 
        intervalFile = fopen(interval_filename.c_str(), "w");

        for (auto e: out_intervals) {
            seg_q_start = block_start + e.query_start;

            if(seg_q_start < curr_q_chr_start || seg_q_start >= curr_q_chr_end){
                q_index = std::upper_bound(chr_start.cbegin(), chr_start.cend(), seg_q_start) - chr_start.cbegin() - 1;
                curr_q_chr_index = q_index;
                curr_q_chr       = chr_name[curr_q_chr_index];
                curr_q_chr_start = chr_start[curr_q_chr_index];
                curr_q_chr_end   = curr_q_chr_start + chr_len[curr_q_chr_index];
            }

            fprintf(intervalFile, "%s\t%lu\t%lu\n", curr_q_chr.c_str(), seg_q_start-curr_q_chr_start, seg_q_start+e.len+1-curr_q_chr_start);
        }

        if(cfg.markend)
            fprintf(intervalFile, "# segalign_repeat_masker end-of-file\n");

        fclose(intervalFile);
    }

    std::get<0>(op).try_put(token);
};
