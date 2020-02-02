#include <atomic>
#include "ntcoding.h"
#include "graph.h"
#include <stddef.h>

std::mutex io_lock;

void segment_printer_body::operator()(printer_input input, printer_node::output_ports_type & op){
    auto &payload = get<0>(input); 
    size_t token = get<1>(input);

    auto &index = get<0>(payload);
    auto &fw_segments = get<1>(payload);
    auto &rc_segments = get<2>(payload);

    std::string filename = "tmp"+std::to_string(index)+".segments";
    std::string maf_filename = "tmp"+std::to_string(index)+".maf";
    std::string filename1 = "tmp"+std::to_string(index)+".segments1";
    FILE* segmentFile = fopen(filename.c_str(), "w");
    int print_header = 1;

    for (auto e: fw_segments) {

        if (print_header == 1) {
            io_lock.lock();
            fprintf(segmentFile, "#name1\tstart1\tend1\tname2\tstart2\tend2\tstrand2\tscore\n");
            print_header = 0;
            io_lock.unlock();
        }

        io_lock.lock();
        fprintf(segmentFile, "ce11.chr1\t%d\t%d\tcb4.chr1\t%d\t%d\t+\t%d\n", (e.ref_start+1), (e.ref_start+e.len), (e.query_start+1), (e.query_start+e.len), e.score);
        io_lock.unlock();
    }

    for (auto e: rc_segments) {

        io_lock.lock();
        fprintf(segmentFile, "ce11.chr1\t%d\t%d\tcb4.chr1\t%d\t%d\t-\t%d\n", (e.ref_start+1), (e.ref_start+e.len), (e.query_start+1), (e.query_start+e.len), e.score);
        io_lock.unlock();
    }

    fclose(segmentFile);

    std::string cmd;
    int status;

    cmd = cfg.lastz_path+" ~/WGA_GPU/data/ce11_chr1.fa ~/WGA_GPU/data/cb4_chr1.fa --format=maf --segments="+filename+" > "+maf_filename;
    status = system(cmd.c_str());
//    cmd = "rm "+filename;
//    status = system(cmd.c_str());

//    cmd = cfg.lastz_path+" ~/WGA_GPU/data/ce11_chr1.fa ~/WGA_GPU/data/cb4_chr1.fa --format=maf --segments="+filename1+" > "+maf_filename;
//    status = system(cmd.c_str());
//    cmd = "rm "+filename1+" "+ filename;
//    status = system(cmd.c_str());

    get<0>(op).try_put(token);
};
