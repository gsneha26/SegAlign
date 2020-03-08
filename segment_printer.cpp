#include "graph.h"

std::mutex io_lock;
int batch = 0;
int batch_old = 0;

void segment_printer_body::operator()(printer_input input, printer_node::output_ports_type & op){
    auto &payload = get<0>(input); 
    size_t token = get<1>(input);

    auto &index = get<0>(payload);
    auto &fw_segments = get<1>(payload);
    auto &rc_segments = get<2>(payload);
    auto &query_chr = get<3>(payload);
    auto &ref_chr = get<4>(payload);

    std::string segment_filename       = "tmp"+std::to_string(index)+"."+ref_chr+"."+query_chr+".segments";
    std::string output_filename   = "tmp"+std::to_string(index)+"."+ref_chr+"."+query_chr+"."+cfg.output_format;

    FILE* segmentFile = fopen(segment_filename.c_str(), "w");

    fprintf(segmentFile, "#name1\tstart1\tend1\tname2\tstart2\tend2\tstrand2\tscore\n");

    for (auto e: fw_segments) {
        fprintf(segmentFile, "%s\t%d\t%d\t%s\t%d\t%d\t+\t%d\n", ref_chr.c_str(),(e.ref_start+1), (e.ref_start+e.len+1), query_chr.c_str(), (e.query_start+1), (e.query_start+e.len+1), e.score);
    }

    for (auto e: rc_segments) {
        fprintf(segmentFile, "%s\t%d\t%d\t%s\t%d\t%d\t-\t%d\n", ref_chr.c_str(), (e.ref_start+1), (e.ref_start+e.len+1), query_chr.c_str(), (e.query_start+1), (e.query_start+e.len+1), e.score);
    }

    fclose(segmentFile);

    std::string cmd;
    int status;

    if(cfg.gapped){
            cmd = "lastz "+cfg.data_folder+"ref/"+ref_chr+".2bit "+cfg.data_folder+"query/"+query_chr+".2bit --format="+ cfg.output_format +" --ydrop="+std::to_string(cfg.ydrop)+" --gappedthresh="+std::to_string(cfg.gappedthresh);
            if(cfg.scoring_file != "")
                cmd = cmd+" --scoring=" + cfg.scoring_file + " --segments="+segment_filename+" > "+output_filename;
            cmd = cmd+" --segments="+segment_filename+" > "+output_filename;
            status = system(cmd.c_str());
    }

    get<0>(op).try_put(token);
};
