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
    uint32_t rc_q_inter_start = q_block_len - cfg.seed_size - q_inter_end;
    uint32_t rc_q_inter_end = q_block_len - cfg.seed_size - q_inter_start;

    std::string base_filename;
    std::string segment_filename;
    std::string output_filename;
    std::string err_filename;
    std::string cmd;

    uint32_t num_hsps = fw_segments.size() + rc_segments.size();

    if(num_hsps > 0){

        size_t start_r_chr = std::upper_bound(r_chr_start.cbegin(), r_chr_start.cend(), r_block_start) - r_chr_start.cbegin() - 1; 
        size_t end_r_chr = std::upper_bound(r_chr_start.cbegin(), r_chr_start.cend(), r_block_end) - r_chr_start.cbegin() - 1; 
        size_t start_q_chr = std::upper_bound(q_chr_start.cbegin(), q_chr_start.cend(), q_block_start+q_inter_start) - q_chr_start.cbegin() - 1; 
        size_t end_q_chr = std::upper_bound(q_chr_start.cbegin(), q_chr_start.cend(), q_block_start+q_inter_end) - q_chr_start.cbegin() - 1; 

        uint32_t num_ref = end_r_chr - start_r_chr + 1;
        uint32_t num_query = end_q_chr - start_q_chr + 1;

        std::vector<std::string> base_files; 

        for(size_t r = start_r_chr; r <= end_r_chr; r++){
            base_files.push_back("tmp"+std::to_string(index)+".block"+std::to_string(block_index)+"."+r_chr_name[r]);
        }

        uint32_t curr_q_chr_index = start_q_chr;
        std::string curr_q_chr = q_chr_name[start_q_chr];
        uint32_t curr_q_chr_file_name = q_chr_file_name[start_q_chr];
        size_t curr_q_chr_start = q_chr_start[start_q_chr]; 
        size_t curr_q_chr_end = curr_q_chr_start + q_chr_len[start_q_chr];

        std::vector <std::vector <std::string>> sorted_segments(num_ref);
        std::string out_str;

        for (auto e: fw_segments) {
            size_t seg_r_start = e.ref_start + r_block_start;
            size_t seg_q_start = q_block_start + e.query_start;
            size_t r_index = std::upper_bound(r_chr_start.cbegin(), r_chr_start.cend(), seg_r_start) - r_chr_start.cbegin() - 1 - start_r_chr;
            size_t q_index = std::upper_bound(q_chr_start.cbegin(), q_chr_start.cend(), seg_q_start) - q_chr_start.cbegin() - 1;

            if(seg_q_start >= curr_q_chr_start && seg_q_start < curr_q_chr_end){
                sorted_segments[r_index].emplace_back(r_chr_name[r_index+start_r_chr] + '\t' + std::to_string(seg_r_start+1-r_chr_start[r_index+start_r_chr]) + '\t' + std::to_string(seg_r_start+e.len+1-r_chr_start[r_index+start_r_chr]) + '\t' + curr_q_chr + '\t' +  std::to_string(seg_q_start+1-curr_q_chr_start) + '\t' + std::to_string(seg_q_start+e.len+1-curr_q_chr_start) + "\t+\t" + std::to_string(e.score) + "\n");
            }
            else{
                for(int i = 0; i < num_ref;i++){
                    if(sorted_segments[i].size() > 0){
                        base_filename = base_files[i]+"."+curr_q_chr+".plus";
                        segment_filename = base_filename+".segments";
                        err_filename = base_filename+".err";
                        output_filename  = base_filename+"."+cfg.output_format;

                        FILE* segmentFile = fopen(segment_filename.c_str(), "w");
                        for(int j = 0; j< sorted_segments[i].size(); j++){
                            fprintf(segmentFile, "%s", sorted_segments[i][j].c_str());
                        }
                        fclose(segmentFile);
                        if(cfg.gapped){

                            cmd = "lastz "+cfg.data_folder+"ref/chr"+std::to_string(r_chr_file_name[i+start_r_chr])+".2bit[nameparse=darkspace] "+cfg.data_folder+"query/chr"+std::to_string(curr_q_chr_file_name)+".2bit[nameparse=darkspace] --format="+ cfg.output_format +" --ydrop="+std::to_string(cfg.ydrop)+" --gappedthresh="+std::to_string(cfg.gappedthresh)+" 2> "+err_filename;
                            if(cfg.notrivial)
                                cmd = cmd+" --notrivial";
                            if(cfg.scoring_file != "")
                                cmd = cmd+" --scoring=" + cfg.scoring_file;
                            cmd = cmd+" --segments="+segment_filename+" --output="+output_filename;

                            io_lock.lock();
                            printf("%s\n", cmd.c_str());
                            io_lock.unlock();
                        }
                        sorted_segments[i].clear();
                    }
                }

                size_t q_index = std::upper_bound(q_chr_start.cbegin(), q_chr_start.cend(), seg_q_start) - q_chr_start.cbegin() - 1;
                curr_q_chr_index = q_index;
                curr_q_chr = q_chr_name[curr_q_chr_index];
                curr_q_chr_file_name = q_chr_file_name[curr_q_chr_index];
                curr_q_chr_start = q_chr_start[curr_q_chr_index];
                curr_q_chr_end = curr_q_chr_start + q_chr_len[curr_q_chr_index];
                sorted_segments[r_index].emplace_back(r_chr_name[r_index+start_r_chr] + '\t' + std::to_string(seg_r_start+1-r_chr_start[r_index+start_r_chr]) + '\t' + std::to_string(seg_r_start+e.len+1-r_chr_start[r_index+start_r_chr]) + '\t' + curr_q_chr + '\t' +  std::to_string(seg_q_start+1-curr_q_chr_start) + '\t' + std::to_string(seg_q_start+e.len+1-curr_q_chr_start) + "\t+\t" + std::to_string(e.score) + "\n");
            }
        }

        for(int i = 0; i < num_ref;i++){
            if(sorted_segments[i].size() > 0){
                base_filename = base_files[i]+"."+curr_q_chr+".plus";
                segment_filename = base_filename+".segments";
                err_filename = base_filename+".err";
                output_filename  = base_filename+"."+cfg.output_format;

                FILE* segmentFile = fopen(segment_filename.c_str(), "w");
                for(int j = 0; j< sorted_segments[i].size(); j++){
                    fprintf(segmentFile, "%s", sorted_segments[i][j].c_str());
                }
                fclose(segmentFile);

                if(cfg.gapped){

                    cmd = "lastz "+cfg.data_folder+"ref/chr"+std::to_string(r_chr_file_name[i+start_r_chr])+".2bit[nameparse=darkspace] "+cfg.data_folder+"query/chr"+std::to_string(curr_q_chr_file_name)+".2bit[nameparse=darkspace] --format="+ cfg.output_format +" --ydrop="+std::to_string(cfg.ydrop)+" --gappedthresh="+std::to_string(cfg.gappedthresh)+" 2> "+err_filename;
                    if(cfg.notrivial)
                        cmd = cmd+" --notrivial";
                    if(cfg.scoring_file != "")
                        cmd = cmd+" --scoring=" + cfg.scoring_file;
                    cmd = cmd+" --segments="+segment_filename+" --output="+output_filename;

                    io_lock.lock();
                    printf("%s\n", cmd.c_str());
                    io_lock.unlock();
                }

                sorted_segments[i].clear();
            }
        }

        start_q_chr = std::upper_bound(rc_q_chr_start.cbegin(), rc_q_chr_start.cend(), q_block_start+rc_q_inter_start) - rc_q_chr_start.cbegin() - 1; 
        end_q_chr = std::upper_bound(rc_q_chr_start.cbegin(), rc_q_chr_start.cend(), q_block_start+rc_q_inter_end) - rc_q_chr_start.cbegin() - 1; 

        curr_q_chr_index = start_q_chr;
        curr_q_chr = rc_q_chr_name[curr_q_chr_index];
        curr_q_chr_file_name = rc_q_chr_file_name[curr_q_chr_index];
        curr_q_chr_start = rc_q_chr_start[curr_q_chr_index];
        curr_q_chr_end = curr_q_chr_start + rc_q_chr_len[curr_q_chr_index];

        
        for (auto e: rc_segments) {
            size_t seg_r_start = e.ref_start + r_block_start;
            size_t seg_q_start = e.query_start + q_block_start;
            size_t r_index = std::upper_bound(r_chr_start.cbegin(), r_chr_start.cend(), seg_r_start) - r_chr_start.cbegin() - 1 - start_r_chr;

            if(seg_q_start >= curr_q_chr_start && seg_q_start < curr_q_chr_end){
                sorted_segments[r_index].emplace_back(r_chr_name[r_index+start_r_chr] + '\t' + std::to_string(seg_r_start+1-r_chr_start[r_index+start_r_chr]) + '\t' + std::to_string(seg_r_start+e.len+1-r_chr_start[r_index+start_r_chr]) + '\t' + curr_q_chr + '\t' +  std::to_string(seg_q_start+1-curr_q_chr_start) + '\t' + std::to_string(seg_q_start+e.len+1-curr_q_chr_start) + "\t-\t" + std::to_string(e.score) + "\n");
            }
            else{
                for(int i = 0; i < num_ref;i++){
                    if(sorted_segments[i].size() > 0){
                        base_filename = base_files[i]+"."+curr_q_chr+".minus";
                        segment_filename = base_filename+".segments";
                        err_filename = base_filename+".err";
                        output_filename  = base_filename+"."+cfg.output_format;

                        FILE* segmentFile = fopen(segment_filename.c_str(), "a");
                        for(int j = 0; j< sorted_segments[i].size(); j++){
                            fprintf(segmentFile, "%s", sorted_segments[i][j].c_str());
                        }
                        fclose(segmentFile);

                        if(cfg.gapped){

                            cmd = "lastz "+cfg.data_folder+"ref/chr"+std::to_string(r_chr_file_name[i+start_r_chr])+".2bit[nameparse=darkspace] "+cfg.data_folder+"query/chr"+std::to_string(curr_q_chr_file_name)+".2bit[nameparse=darkspace] --format="+ cfg.output_format +" --ydrop="+std::to_string(cfg.ydrop)+" --gappedthresh="+std::to_string(cfg.gappedthresh)+" 2> "+err_filename;
                            if(cfg.notrivial)
                                cmd = cmd+" --notrivial";
                            if(cfg.scoring_file != "")
                                cmd = cmd+" --scoring=" + cfg.scoring_file;
                            cmd = cmd+" --segments="+segment_filename+" --output="+output_filename;

                            io_lock.lock();
                            printf("%s\n", cmd.c_str());
                            io_lock.unlock();
                        }
                        sorted_segments[i].clear();
                        
                    }
                }

                size_t q_index = std::upper_bound(rc_q_chr_start.cbegin(), rc_q_chr_start.cend(), seg_q_start) - rc_q_chr_start.cbegin() - 1;
                curr_q_chr_index = q_index;
                curr_q_chr = rc_q_chr_name[curr_q_chr_index];
                curr_q_chr_file_name = rc_q_chr_file_name[curr_q_chr_index];
                curr_q_chr_start = rc_q_chr_start[curr_q_chr_index];
                curr_q_chr_end = curr_q_chr_start + rc_q_chr_len[curr_q_chr_index];
                sorted_segments[r_index].emplace_back(r_chr_name[r_index+start_r_chr] + '\t' + std::to_string(seg_r_start+1-r_chr_start[r_index+start_r_chr]) + '\t' + std::to_string(seg_r_start+e.len+1-r_chr_start[r_index+start_r_chr]) + '\t' + curr_q_chr + '\t' +  std::to_string(seg_q_start+1-curr_q_chr_start) + '\t' + std::to_string(seg_q_start+e.len+1-curr_q_chr_start) + "\t-\t" + std::to_string(e.score) + "\n");
            }
        }

        for(int i = 0; i < num_ref;i++){
            if(sorted_segments[i].size() > 0){
                base_filename = base_files[i]+"."+curr_q_chr+".minus";
                segment_filename = base_filename+".segments";
                err_filename = base_filename+".err";
                output_filename  = base_filename+"."+cfg.output_format;

                FILE* segmentFile = fopen(segment_filename.c_str(), "a");
                for(int j = 0; j< sorted_segments[i].size(); j++){
                    fprintf(segmentFile, "%s", sorted_segments[i][j].c_str());
                }
                fclose(segmentFile);
                if(cfg.gapped){

                    cmd = "lastz "+cfg.data_folder+"ref/chr"+std::to_string(r_chr_file_name[i+start_r_chr])+".2bit[nameparse=darkspace] "+cfg.data_folder+"query/chr"+std::to_string(curr_q_chr_file_name)+".2bit[nameparse=darkspace] --format="+ cfg.output_format +" --ydrop="+std::to_string(cfg.ydrop)+" --gappedthresh="+std::to_string(cfg.gappedthresh)+" 2> "+err_filename;
                    if(cfg.notrivial)
                        cmd = cmd+" --notrivial";
                    if(cfg.scoring_file != "")
                        cmd = cmd+" --scoring=" + cfg.scoring_file;
                    cmd = cmd+" --segments="+segment_filename+" --output="+output_filename;

                    io_lock.lock();
                    printf("%s\n", cmd.c_str());
                    io_lock.unlock();
                }
                sorted_segments[i].clear();
            }
        }
    }

    get<0>(op).try_put(token);
};
