#include <algorithm>
#include <iterator>
#include <string>
#include <map>
#include <mutex>
#include "graph.h"
#include "store.h"

std::mutex io_lock;

void segment_printer_body::operator()(printer_input input, printer_node::output_ports_type & op){

    auto &payload = std::get<0>(input);
    size_t token  = std::get<1>(input);

    auto &seeder_payload = std::get<0>(payload);
    auto &fw_hsps    = std::get<1>(payload);
    auto &rc_hsps    = std::get<2>(payload);

    auto &block_data    = std::get<0>(seeder_payload);
    auto &interval_data = std::get<1>(seeder_payload);

    int r_block_index = block_data.r_index - 1;
    int q_block_index = block_data.q_index;
    size_t r_block_start = block_data.r_start;
    size_t q_block_start = block_data.q_start;
    uint32_t r_block_len = block_data.r_len;
    uint32_t q_block_len = block_data.q_len;

    uint32_t q_inter_start = interval_data.start;
    uint32_t q_inter_end   = interval_data.end;
    uint32_t num_invoked   = interval_data.num_invoked;

    uint32_t index = num_invoked;
    size_t r_block_end = r_block_start + r_block_len - 1;
    size_t rc_q_inter_start = q_block_len - q_inter_end;
    size_t rc_q_inter_end = q_block_len - q_inter_start;

    std::string base_filename;
    std::string segment_filename;
    std::string output_filename;
    std::string err_filename;
    std::string cmd;

    uint32_t fw_num_hsps  = fw_hsps.size();
    uint32_t rc_num_hsps  = rc_hsps.size();
    uint32_t num_hsps = fw_num_hsps + rc_num_hsps;

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

        size_t seg_r_start;
        size_t seg_q_start;
        size_t r_index;
        size_t q_index;

        std::string out_str;
        FILE* segmentFile;

        if(fw_num_hsps > 0){

            base_filename = "tmp"+std::to_string(index)+".block"+std::to_string(q_block_index)+".r"+std::to_string(r_block_start)+".plus"; 
            segment_filename = base_filename+".segments";
            segmentFile = fopen(segment_filename.c_str(), "w");

            for (auto e: fw_hsps) {
                seg_r_start = e.ref_start + r_block_start;
                seg_q_start = q_block_start + e.query_start;
                r_index = std::upper_bound(r_chr_start.cbegin(), r_chr_start.cend(), seg_r_start) - r_chr_start.cbegin() - 1;

                if(seg_q_start < curr_q_chr_start || seg_q_start >= curr_q_chr_end){
                    q_index = std::upper_bound(q_chr_start.cbegin(), q_chr_start.cend(), seg_q_start) - q_chr_start.cbegin() - 1;
                    curr_q_chr_index = q_index;
                    curr_q_chr = q_chr_name[curr_q_chr_index];
                    curr_q_chr_file_name = q_chr_file_name[curr_q_chr_index];
                    curr_q_chr_start = q_chr_start[curr_q_chr_index];
                    curr_q_chr_end = curr_q_chr_start + q_chr_len[curr_q_chr_index];
                }

                out_str = r_chr_name[r_index] + '\t' + std::to_string(seg_r_start+1-r_chr_start[r_index]) + '\t' + std::to_string(seg_r_start+e.len+1-r_chr_start[r_index]) + '\t' + curr_q_chr + '\t' +  std::to_string(seg_q_start+1-curr_q_chr_start) + '\t' + std::to_string(seg_q_start+e.len+1-curr_q_chr_start) + "\t+\t" + std::to_string(e.score) + "\n";
                fprintf(segmentFile, "%s", out_str.c_str());
            }

            fclose(segmentFile);

            if(cfg.gapped){

                err_filename = base_filename+".err";
                output_filename  = base_filename+"."+cfg.output_format;

                cmd = "lastz "+cfg.data_folder+"ref.2bit[nameparse=darkspace][multiple][subset=ref_block"+std::to_string(r_block_index)+".name] "+cfg.data_folder+"query.2bit[nameparse=darkspace][subset=query_block"+std::to_string(q_block_index)+".name] --format="+ cfg.output_format +" --ydrop="+std::to_string(cfg.ydrop)+" --gappedthresh="+std::to_string(cfg.gappedthresh)+" --strand=plus";
                if(cfg.ambiguous != "")
                    cmd = cmd+" --ambiguous="+cfg.ambiguous;
                if(cfg.notrivial)
                    cmd = cmd+" --notrivial";
                if(cfg.scoring_file != "")
                    cmd = cmd+" --scores=" + cfg.scoring_file;
                cmd = cmd+" --segments="+segment_filename+" --output="+output_filename+" 2> "+err_filename;

                io_lock.lock();
                printf("%s\n", cmd.c_str());
                io_lock.unlock();
            }
        }

        start_q_chr = std::upper_bound(rc_q_chr_start.cbegin(), rc_q_chr_start.cend(), q_block_start+rc_q_inter_start) - rc_q_chr_start.cbegin() - 1; 
        end_q_chr = std::upper_bound(rc_q_chr_start.cbegin(), rc_q_chr_start.cend(), q_block_start+rc_q_inter_end) - rc_q_chr_start.cbegin() - 1; 

        curr_q_chr_index = end_q_chr;
        curr_q_chr = rc_q_chr_name[curr_q_chr_index];
        curr_q_chr_file_name = rc_q_chr_file_name[curr_q_chr_index];
        curr_q_chr_start = rc_q_chr_start[curr_q_chr_index];
        curr_q_chr_end = curr_q_chr_start + rc_q_chr_len[curr_q_chr_index];

        if(rc_num_hsps > 0){
            base_filename = "tmp"+std::to_string(index)+".block"+std::to_string(q_block_index)+".r"+std::to_string(r_block_start)+".minus"; 
            segment_filename = base_filename+".segments";
            segmentFile = fopen(segment_filename.c_str(), "w");

            for(int r = rc_hsps.size()-1; r >= 0; r--){
                auto e =  rc_hsps[r];
                seg_r_start = e.ref_start + r_block_start;
                seg_q_start = e.query_start + q_block_start;
                r_index = std::upper_bound(r_chr_start.cbegin(), r_chr_start.cend(), seg_r_start) - r_chr_start.cbegin() - 1;

                if(seg_q_start < curr_q_chr_start || seg_q_start >= curr_q_chr_end){
                    q_index = std::upper_bound(rc_q_chr_start.cbegin(), rc_q_chr_start.cend(), seg_q_start) - rc_q_chr_start.cbegin() - 1;
                    curr_q_chr_index = q_index;
                    curr_q_chr = rc_q_chr_name[curr_q_chr_index];
                    curr_q_chr_file_name = rc_q_chr_file_name[curr_q_chr_index];
                    curr_q_chr_start = rc_q_chr_start[curr_q_chr_index];
                    curr_q_chr_end = curr_q_chr_start + rc_q_chr_len[curr_q_chr_index];
                }

                out_str = r_chr_name[r_index] + '\t' + std::to_string(seg_r_start+1-r_chr_start[r_index]) + '\t' + std::to_string(seg_r_start+e.len+1-r_chr_start[r_index]) + '\t' + curr_q_chr + '\t' +  std::to_string(seg_q_start+1-curr_q_chr_start) + '\t' + std::to_string(seg_q_start+e.len+1-curr_q_chr_start) + "\t-\t" + std::to_string(e.score) + "\n";
                fprintf(segmentFile, "%s", out_str.c_str());
            }

            fclose(segmentFile);

            if(cfg.gapped){

                err_filename = base_filename+".err";
                output_filename  = base_filename+"."+cfg.output_format;

                cmd = "lastz "+cfg.data_folder+"ref.2bit[nameparse=darkspace][multiple][subset=ref_block"+std::to_string(r_block_index)+".name] "+cfg.data_folder+"query.2bit[nameparse=darkspace][subset=query_block"+std::to_string(q_block_index)+".name] --format="+ cfg.output_format +" --ydrop="+std::to_string(cfg.ydrop)+" --gappedthresh="+std::to_string(cfg.gappedthresh)+" --strand=minus";
                if(cfg.ambiguous != "")
                    cmd = cmd+" --ambiguous="+cfg.ambiguous;
                if(cfg.notrivial)
                    cmd = cmd+" --notrivial";
                if(cfg.scoring_file != "")
                    cmd = cmd+" --scores=" + cfg.scoring_file;
                cmd = cmd+" --segments="+segment_filename+" --output="+output_filename+" 2> "+err_filename;

                io_lock.lock();
                printf("%s\n", cmd.c_str());
                io_lock.unlock();
            }
        }
    }

    std::get<0>(op).try_put(token);
};
