#include "graph.h"
#include <algorithm>
#include <iterator>
#include <iostream>
#include <map>

std::mutex io_lock;

void segment_printer_body::operator()(printer_input input, printer_node::output_ports_type & op){

    auto &payload = get<0>(input); 
    size_t token  = get<1>(input);

    auto &index = get<0>(payload);
    auto &fw_segments = get<1>(payload);
    auto &rc_segments = get<2>(payload);
    auto &block_index = get<3>(payload);
    auto &r_block_start = get<4>(payload);
    auto &r_block_end = get<5>(payload);
    auto &q_block_start = get<6>(payload);
    auto &q_inter_start = get<7>(payload);
    auto &q_inter_end = get<8>(payload);
    auto &q_block_len = get<9>(payload);
    auto &r_block_index = get<10>(payload);
    r_block_index--;
    uint32_t rc_q_inter_start = q_block_len - cfg.seed_size - q_inter_end;
    uint32_t rc_q_inter_end = q_block_len - cfg.seed_size - q_inter_start;

    std::string base_filename;
    std::string segment_filename;
    std::string output_filename;
    std::string err_filename;
    std::string cmd;

    uint32_t num_fw_hsps = fw_segments.size();
    uint32_t num_rc_hsps = rc_segments.size();
    uint32_t num_hsps = num_fw_hsps + num_rc_hsps;

    if(num_hsps > 0){

        size_t start_r_chr = std::upper_bound(r_chr_start.cbegin(), r_chr_start.cend(), r_block_start) - r_chr_start.cbegin() - 1; 
        size_t end_r_chr = std::upper_bound(r_chr_start.cbegin(), r_chr_start.cend(), r_block_end) - r_chr_start.cbegin() - 1; 
        size_t start_q_chr = std::upper_bound(q_chr_start.cbegin(), q_chr_start.cend(), q_block_start+q_inter_start) - q_chr_start.cbegin() - 1; 
        size_t end_q_chr = std::upper_bound(q_chr_start.cbegin(), q_chr_start.cend(), q_block_start+q_inter_end) - q_chr_start.cbegin() - 1; 

        uint32_t curr_q_chr_index = start_q_chr;
        std::string curr_q_chr = q_chr_name[start_q_chr];
        uint32_t curr_q_chr_file_name = q_chr_file_name[start_q_chr];
        size_t curr_q_chr_start = q_chr_start[start_q_chr]; 
        size_t curr_q_chr_end = curr_q_chr_start + q_chr_len[start_q_chr];

        std::string out_str;
        FILE* segmentFile;

        uint32_t file_index = 0;
        uint32_t hsps_printed = 0;

        if(num_fw_hsps > 0){

            base_filename = "tmp"+std::to_string(index)+".file"+std::to_string(file_index)+".block"+std::to_string(block_index)+".r"+std::to_string(r_block_start)+".plus"; 
            segment_filename = base_filename+".segments";
            segmentFile = fopen(segment_filename.c_str(), "w");

            for (auto e: fw_segments) {
                if(hsps_printed > MAX_HSPS_FILE){
                    fclose(segmentFile);

                    cmd = "aws s3 cp "+segment_filename+" "+cfg.data_folder;

                    io_lock.lock();
                    printf("%s\n", cmd.c_str());
                    io_lock.unlock();

                    file_index += 1;
                    hsps_printed = 1;
                    base_filename = "tmp"+std::to_string(index)+".file"+std::to_string(file_index)+".block"+std::to_string(block_index)+".r"+std::to_string(r_block_start)+".plus"; 
                    segment_filename = base_filename+".segments";
                    segmentFile = fopen(segment_filename.c_str(), "w");
                }
                else{
                    hsps_printed += 1;
                }

                size_t seg_r_start = e.ref_start + r_block_start;
                size_t seg_q_start = q_block_start + e.query_start;
                size_t r_index = std::upper_bound(r_chr_start.cbegin(), r_chr_start.cend(), seg_r_start) - r_chr_start.cbegin() - 1;

                if(seg_q_start >= curr_q_chr_start && seg_q_start < curr_q_chr_end){
                    out_str = r_chr_name[r_index] + '\t' + std::to_string(seg_r_start+1-r_chr_start[r_index]) + '\t' + std::to_string(seg_r_start+e.len+1-r_chr_start[r_index]) + '\t' + curr_q_chr + '\t' +  std::to_string(seg_q_start+1-curr_q_chr_start) + '\t' + std::to_string(seg_q_start+e.len+1-curr_q_chr_start) + "\t+\t" + std::to_string(e.score) + "\n";
                }
                else{
                    size_t q_index = std::upper_bound(q_chr_start.cbegin(), q_chr_start.cend(), seg_q_start) - q_chr_start.cbegin() - 1;
                    curr_q_chr_index = q_index;
                    curr_q_chr = q_chr_name[curr_q_chr_index];
                    curr_q_chr_file_name = q_chr_file_name[curr_q_chr_index];
                    curr_q_chr_start = q_chr_start[curr_q_chr_index];
                    curr_q_chr_end = curr_q_chr_start + q_chr_len[curr_q_chr_index];

                    out_str = r_chr_name[r_index] + '\t' + std::to_string(seg_r_start+1-r_chr_start[r_index]) + '\t' + std::to_string(seg_r_start+e.len+1-r_chr_start[r_index]) + '\t' + curr_q_chr + '\t' +  std::to_string(seg_q_start+1-curr_q_chr_start) + '\t' + std::to_string(seg_q_start+e.len+1-curr_q_chr_start) + "\t+\t" + std::to_string(e.score) + "\n";
                }
                fprintf(segmentFile, "%s", out_str.c_str());
            }

            fclose(segmentFile);

            cmd = "aws s3 cp "+segment_filename+" "+cfg.data_folder;

            io_lock.lock();
            printf("%s\n", cmd.c_str());
            io_lock.unlock();

        }

        start_q_chr = std::upper_bound(rc_q_chr_start.cbegin(), rc_q_chr_start.cend(), q_block_start+rc_q_inter_start) - rc_q_chr_start.cbegin() - 1; 
        end_q_chr = std::upper_bound(rc_q_chr_start.cbegin(), rc_q_chr_start.cend(), q_block_start+rc_q_inter_end) - rc_q_chr_start.cbegin() - 1; 

        curr_q_chr_index = end_q_chr;
        curr_q_chr = rc_q_chr_name[curr_q_chr_index];
        curr_q_chr_file_name = rc_q_chr_file_name[curr_q_chr_index];
        curr_q_chr_start = rc_q_chr_start[curr_q_chr_index];
        curr_q_chr_end = curr_q_chr_start + rc_q_chr_len[curr_q_chr_index];

        file_index = 0;
        hsps_printed = 0;

        if(num_rc_hsps > 0){
            base_filename = "tmp"+std::to_string(index)+".file"+std::to_string(file_index)+".block"+std::to_string(block_index)+".r"+std::to_string(r_block_start)+".minus"; 
            segment_filename = base_filename+".segments";
            segmentFile = fopen(segment_filename.c_str(), "w");

            for(int r = num_rc_hsps-1; r >= 0; r--){
                if(hsps_printed > MAX_HSPS_FILE){
                    fclose(segmentFile);

                    cmd = "aws s3 cp "+segment_filename+" "+cfg.data_folder;

                    io_lock.lock();
                    printf("%s\n", cmd.c_str());
                    io_lock.unlock();

                    file_index += 1;
                    hsps_printed = 1;
                    base_filename = "tmp"+std::to_string(index)+".file"+std::to_string(file_index)+".block"+std::to_string(block_index)+".r"+std::to_string(r_block_start)+".minus"; 
                    segment_filename = base_filename+".segments";
                    segmentFile = fopen(segment_filename.c_str(), "w");
                }
                else{
                    hsps_printed += 1;
                }
                auto e =  rc_segments[r];
                size_t seg_r_start = e.ref_start + r_block_start;
                size_t seg_q_start = e.query_start + q_block_start;
                size_t r_index = std::upper_bound(r_chr_start.cbegin(), r_chr_start.cend(), seg_r_start) - r_chr_start.cbegin() - 1;

                if(seg_q_start >= curr_q_chr_start && seg_q_start < curr_q_chr_end){
                    out_str = r_chr_name[r_index] + '\t' + std::to_string(seg_r_start+1-r_chr_start[r_index]) + '\t' + std::to_string(seg_r_start+e.len+1-r_chr_start[r_index]) + '\t' + curr_q_chr + '\t' +  std::to_string(seg_q_start+1-curr_q_chr_start) + '\t' + std::to_string(seg_q_start+e.len+1-curr_q_chr_start) + "\t-\t" + std::to_string(e.score) + "\n";
                }
                else{

                    size_t q_index = std::upper_bound(rc_q_chr_start.cbegin(), rc_q_chr_start.cend(), seg_q_start) - rc_q_chr_start.cbegin() - 1;
                    curr_q_chr_index = q_index;
                    curr_q_chr = rc_q_chr_name[curr_q_chr_index];
                    curr_q_chr_file_name = rc_q_chr_file_name[curr_q_chr_index];
                    curr_q_chr_start = rc_q_chr_start[curr_q_chr_index];
                    curr_q_chr_end = curr_q_chr_start + rc_q_chr_len[curr_q_chr_index];

                    out_str = r_chr_name[r_index] + '\t' + std::to_string(seg_r_start+1-r_chr_start[r_index]) + '\t' + std::to_string(seg_r_start+e.len+1-r_chr_start[r_index]) + '\t' + curr_q_chr + '\t' +  std::to_string(seg_q_start+1-curr_q_chr_start) + '\t' + std::to_string(seg_q_start+e.len+1-curr_q_chr_start) + "\t-\t" + std::to_string(e.score) + "\n";
                }
                fprintf(segmentFile, "%s", out_str.c_str());
            }

            fclose(segmentFile);

            cmd = "aws s3 cp "+segment_filename+" "+cfg.data_folder;

            io_lock.lock();
            printf("%s\n", cmd.c_str());
            io_lock.unlock();

        }
    }

    get<0>(op).try_put(token);
};
