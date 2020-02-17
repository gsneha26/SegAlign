#include <atomic>
#include "ntcoding.h"
#include "graph.h"
#include <stddef.h>

std::mutex io_lock;
//int* diagHashminus = (int*) calloc(1000000000, sizeof(int));
//int* diagHashplus = (int*) calloc(1000000000, sizeof(int));
int batch = 0;
int batch_old = 0;

void segment_printer_body::operator()(printer_input input, printer_node::output_ports_type & op){
    auto &payload = get<0>(input); 
    size_t token = get<1>(input);

    auto &index = get<0>(payload);
    auto &fw_segments = get<1>(payload);
    auto &rc_segments = get<2>(payload);
    auto &query_chr = get<3>(payload);

    std::string filename       = "tmp"+std::to_string(index)+"."+query_chr+".segments";
    std::string maf_filename   = "tmp"+std::to_string(index)+"."+query_chr+".maf";

    FILE* segmentFile = fopen(filename.c_str(), "w");

    fprintf(segmentFile, "#name1\tstart1\tend1\tname2\tstart2\tend2\tstrand2\tscore\n");

    for (auto e: fw_segments) {
        fprintf(segmentFile, "chr1\t%d\t%d\t%s\t%d\t%d\t+\t%d\n", (e.ref_start+1), (e.ref_start+e.len+1), query_chr.c_str(), (e.query_start+1), (e.query_start+e.len+1), e.score);
    }

    for (auto e: rc_segments) {
        fprintf(segmentFile, "chr1\t%d\t%d\t%s\t%d\t%d\t-\t%d\n", (e.ref_start+1), (e.ref_start+e.len+1), query_chr.c_str(), (e.query_start+1), (e.query_start+e.len+1), e.score);
    }

//        int diag = 0;
//        int diag_old = 0;
//        int r_extent = 0;
//        int q_extent = 0;
//        
//        for (auto e: fw_segments) {
//        
//            diag = e.ref_start-e.query_start;
//            if(diag != diag_old){
//                r_extent = 0;
//                q_extent = 0;
//        
//                if(e.ref_start > r_extent || e.query_start > q_extent){
//                    fprintf(segmentFile, "%s.chr1\t%d\t%d\t%s.chr1\t%d\t%d\t+\t%d\n", cfg.reference_name.c_str(), (e.ref_start+1), (e.ref_start+e.len+1), cfg.query_name.c_str(), (e.query_start+1), (e.query_start+e.len+1), e.score);
//                    r_extent = e.ref_start + e.len;
//                    q_extent = e.query_start + e.len;
//                    diag_old = diag;
//                }
//            }
//        }
//
//        diag = 0;
//        diag_old = 0;
//        r_extent = 0;
//        q_extent = 0;
//        
//        for (auto e: rc_segments) {
//        
//            diag = e.ref_start-e.query_start;
//            if(diag != diag_old){
//                r_extent = 0;
//                q_extent = 0;
//        
//                if(e.ref_start > r_extent || e.query_start > q_extent){
//                    fprintf(segmentFile, "%s.chr1\t%d\t%d\t%s.chr1\t%d\t%d\t-\t%d\n", cfg.reference_name.c_str(), (e.ref_start+1), (e.ref_start+e.len+1), cfg.query_name.c_str(), (e.query_start+1), (e.query_start+e.len+1), e.score);
//                    r_extent = e.ref_start + e.len;
//                    q_extent = e.query_start + e.len;
//                    diag_old = diag;
//                }
//            }
//        }

//        int diag;
//        for (auto e: fw_segments) {
//            diag = e.ref_start-e.query_start+cfg.query_len;
//            if(e.query_start > diagHashplus[diag]){
//                fprintf(segmentFile, "%s.chr1\t%d\t%d\t%s.chr1\t%d\t%d\t+\t%d\n", cfg.reference_name.c_str(), (e.ref_start+1), (e.ref_start+e.len+1), cfg.query_name.c_str(), (e.query_start+1), (e.query_start+e.len+1), e.score);
//                diagHashplus[diag] = e.query_start+e.len;
//            }
//        }
//
//        diag = 0;
//        for (auto e: rc_segments) {
//            diag = e.ref_start-e.query_start+cfg.query_len;
//            if(e.query_start > diagHashminus[diag]){
//                fprintf(segmentFile, "%s.chr1\t%d\t%d\t%s.chr1\t%d\t%d\t-\t%d\n", cfg.reference_name.c_str(), (e.ref_start+1), (e.ref_start+e.len+1), cfg.query_name.c_str(), (e.query_start+1), (e.query_start+e.len+1), e.score);
//                diagHashminus[diag] = e.query_start+e.len;
//            }
//        }

    fclose(segmentFile);

    std::string cmd;
    int status;

    if(cfg.do_gapped){
        cmd = cfg.lastz_path+" "+cfg.data_folder+cfg.reference_name+"/chr1.2bit "+cfg.data_folder+cfg.query_name+"/"+query_chr+".2bit --format=maf- --segments="+filename+" > "+maf_filename;
        status = system(cmd.c_str());
    }
//    cmd = "rm "+filename_plus+" "+filename_minus;
//    status = system(cmd.c_str());

    get<0>(op).try_put(token);
};
