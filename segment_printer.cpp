#include "graph.h"
#include <algorithm>
#include <iterator>
#include <vector>
#include <iostream>

std::mutex io_lock;

void segment_printer_body::operator()(printer_input input, printer_node::output_ports_type & op){

    auto &payload = get<0>(input); 
    size_t token  = get<1>(input);

    auto &index = get<0>(payload);
    auto &fw_segments = get<1>(payload);
    auto &rc_segments = get<2>(payload);
    auto &block_index = get<3>(payload);
    auto &r_block_start = get<4>(payload);
    auto &q_block_start = get<5>(payload);

    std::string base_filename;
    std::string segment_filename;
    std::string output_filename;
    std::string err_filename;
    std::string cmd;

    for (auto e: fw_segments) {
        size_t seg_r_start = e.ref_start + r_block_start;
        size_t seg_q_start = e.query_start + q_block_start;
        size_t r_index = std::upper_bound(r_chr_start.cbegin(), r_chr_start.cend(), seg_r_start) - r_chr_start.cbegin() - 1;
        size_t q_index = std::upper_bound(q_chr_start.cbegin(), q_chr_start.cend(), seg_q_start) - q_chr_start.cbegin() - 1;
        std::cout << seg_r_start << " " << r_chr_name[r_index] << " " << seg_r_start-r_chr_start[r_index] +1 << " " << seg_r_start-r_chr_start[r_index] +1+e.len<< " " << seg_q_start << " " << q_chr_name[q_index] << " " << seg_q_start-q_chr_start[q_index] +1<< " " << seg_q_start-q_chr_start[q_index]+1+e.len << " " << e.score << std::endl;
    }

    /*
    base_filename = "tmp"+std::to_string(index);//+".chr"+std::to_string(r_index)+".chr"+std::to_string(q_index);
    segment_filename = base_filename+".segments";
    err_filename = base_filename+".err";
    output_filename  = base_filename+"."+cfg.output_format;

    FILE* segmentFile = fopen(segment_filename.c_str(), "w");

    if(cfg.gapped){
        fprintf(segmentFile, "#name1\tstart1\tend1\tname2\tstart2\tend2\tstrand2\tscore\n");
    }

    for (auto e: fw_segments) {
        fprintf(segmentFile, "%s\t%d\t%d\t%s\t%d\t%d\t+\t%d\n", ref_chr.c_str(),(e.ref_start+1), (e.ref_start+e.len+1), query_chr.c_str(), (e.query_start+1), (e.query_start+e.len+1), e.score);
    }

    for (auto e: rc_segments) {
        fprintf(segmentFile, "%s\t%d\t%d\t%s\t%d\t%d\t-\t%d\n", ref_chr.c_str(), (e.ref_start+1), (e.ref_start+e.len+1), query_chr.c_str(), (e.query_start+1), (e.query_start+e.len+1), e.score);
    }

    fclose(segmentFile);

    if(cfg.gapped){

        std::string cmd;

        cmd = "lastz "+cfg.data_folder+"ref/chr"+std::to_string(index)+".2bit[nameparse=darkspace] "+cfg.data_folder+"query/chr"+std::to_string(index)+".2bit[nameparse=darkspace] --format="+ cfg.output_format +" --ydrop="+std::to_string(cfg.ydrop)+" --gappedthresh="+std::to_string(cfg.gappedthresh)+" 2> "+err_filename;
        if(cfg.notrivial)
            cmd = cmd+" --notrivial";
        if(cfg.scoring_file != "")
            cmd = cmd+" --scoring=" + cfg.scoring_file;
        cmd = cmd+" --segments="+segment_filename+" --output="+output_filename;
            
	io_lock.lock();
	printf("%s\n", cmd.c_str());
	io_lock.unlock();
    }
    */

    get<0>(op).try_put(token);
};
