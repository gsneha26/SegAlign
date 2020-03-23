#include "graph.h"

std::mutex io_lock;

void segment_printer_body::operator()(printer_input input, printer_node::output_ports_type & op){

    auto &payload = get<0>(input); 
    size_t token = get<1>(input);

    auto &index = get<0>(payload);
    auto &fw_segments = get<1>(payload);
    auto &rc_segments = get<2>(payload);
    auto &query_chr = get<3>(payload);
    auto &ref_chr = get<4>(payload);
    auto &r_index = get<5>(payload);
    auto &q_index = get<6>(payload);

    std::string segment_filename;
    std::string output_filename;
    std::string cmd;

    segment_filename = "tmp"+std::to_string(index)+".ref"+std::to_string(r_index)+".query"+std::to_string(q_index)+".segments";
    output_filename  = "tmp"+std::to_string(index)+".ref"+std::to_string(r_index)+".query"+std::to_string(q_index)+"."+cfg.output_format;

    FILE* segmentFile = fopen(segment_filename.c_str(), "w");

    fprintf(segmentFile, "#name1\tstart1\tend1\tname2\tstart2\tend2\tstrand2\tscore\n");

    for (auto e: fw_segments) {
        fprintf(segmentFile, "%s\t%d\t%d\t%s\t%d\t%d\t+\t%d\n", ref_chr.c_str(),(e.ref_start+1), (e.ref_start+e.len+1), query_chr.c_str(), (e.query_start+1), (e.query_start+e.len+1), e.score);
    }

    for (auto e: rc_segments) {
        fprintf(segmentFile, "%s\t%d\t%d\t%s\t%d\t%d\t-\t%d\n", ref_chr.c_str(), (e.ref_start+1), (e.ref_start+e.len+1), query_chr.c_str(), (e.query_start+1), (e.query_start+e.len+1), e.score);
    }

    fclose(segmentFile);

    if(cfg.gapped){

        std::string cmd;

        cmd = "lastz "+cfg.data_folder+"ref/ref"+std::to_string(r_index)+".2bit[nameparse=darkspace] "+cfg.data_folder+"query/query"+std::to_string(q_index)+".2bit[nameparse=darkspace] --format="+ cfg.output_format +" --ydrop="+std::to_string(cfg.ydrop)+" --gappedthresh="+std::to_string(cfg.gappedthresh);
        if(cfg.notrivial)
            cmd = cmd+" --notrivial";
        if(cfg.scoring_file != "")
            cmd = cmd+" --scoring=" + cfg.scoring_file;
        cmd = cmd+" --segments="+segment_filename+" --output="+output_filename;
            
	io_lock.lock();
	printf("%s\n", cmd.c_str());
	io_lock.unlock();
    }

    get<0>(op).try_put(token);
};
