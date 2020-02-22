/*
MIT License

Copyright (c) 2019 Sneha D. Goenka, Yatish Turakhia, Gill Bejerano and William J. Dally

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <assert.h>
#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <tbb/task_scheduler_init.h>

#include <zlib.h>
#include <algorithm>
#include <iterator>
#include <vector>
#include "ConfigFile.h"
#include "graph.h"
#include "kseq.h"
#include "DRAM.h"

////////////////////////////////////////////////////////////////////////////////

KSEQ_INIT2(, gzFile, gzread)

struct timeval start_time, end_time, start_time1;
long useconds, seconds, mseconds;

Configuration cfg;
SeedPosTable *sa;

// query
std::vector<std::string> q_chr_id;
std::vector<uint32_t>  q_chr_len;
std::vector<size_t>  q_chr_coord;
int a = 0;

////////////////////////////////////////////////////////////////////////////////
char* RevComp(bond::blob read) {

    size_t read_len = read.size();

    char* rc = (char*) malloc(read_len*sizeof(char));

    char* seq = (char*)read.data();
    for (size_t r = 0, i = read_len; i-- > 0;) {
        if (seq[i] != 'a' && seq[i] != 'A' &&
                seq[i] != 'c' && seq[i] != 'C' &&
                seq[i] != 'g' && seq[i] != 'G' &&
                seq[i] != 't' && seq[i] != 'T' &&
                seq[i] != 'n' && seq[i] != 'N') {
            fprintf(stderr, "Bad Nt char: %c\n", seq[i]);
            exit(1);
        }
        switch (seq[i]) {
            case 'a': rc[r++] = 't';
                      break;

            case 'A': rc[r++] = 'T';
                      break;

            case 'c': rc[r++] = 'g';
                      break;

            case 'C': rc[r++] = 'G';
                      break;

            case 'g': rc[r++] = 'c';
                      break;

            case 'G': rc[r++] = 'C';
                      break;

            case 't': rc[r++] = 'a';
                      break;

            case 'T': rc[r++] = 'A';
                      break;

            case 'n': rc[r++] = 'n';
                      break;

            case 'N': rc[r++] = 'N';
                      break;
        }
    }

    return rc;
}

int main(int argc, char** argv){

    if (argc != 1) {
        printf("Usage: %s \n", argv[0]);
        return EXIT_FAILURE;
    }

    ConfigFile cfg_file("params.cfg");
    gettimeofday(&start_time, NULL);
    fprintf(stderr, "Loading configuration file ...");

    // FASTA files
    cfg.reference_name      = (std::string) cfg_file.Value("FASTA_files", "reference_name"); 
    cfg.query_name          = (std::string) cfg_file.Value("FASTA_files", "query_name"); 
    cfg.reference_filename  = (std::string) cfg_file.Value("FASTA_files", "reference_filename"); 
    cfg.query_filename      = (std::string) cfg_file.Value("FASTA_files", "query_filename"); 
    cfg.data_folder         = (std::string) cfg_file.Value("FASTA_files", "data_folder");

    // GACT scoring
    cfg.gact_sub_mat[0]     = cfg_file.Value("Scoring", "sub_AA");
    cfg.gact_sub_mat[1]     = cfg_file.Value("Scoring", "sub_AC");
    cfg.gact_sub_mat[2]     = cfg_file.Value("Scoring", "sub_AG");
    cfg.gact_sub_mat[3]     = cfg_file.Value("Scoring", "sub_AT");
    cfg.gact_sub_mat[4]     = cfg_file.Value("Scoring", "sub_CC");
    cfg.gact_sub_mat[5]     = cfg_file.Value("Scoring", "sub_CG");
    cfg.gact_sub_mat[6]     = cfg_file.Value("Scoring", "sub_CT");
    cfg.gact_sub_mat[7]     = cfg_file.Value("Scoring", "sub_GG");
    cfg.gact_sub_mat[8]     = cfg_file.Value("Scoring", "sub_GT");
    cfg.gact_sub_mat[9]     = cfg_file.Value("Scoring", "sub_TT");
    cfg.gact_sub_mat[10]    = cfg_file.Value("Scoring", "sub_N");
    cfg.gap_open            = cfg_file.Value("Scoring", "gap_open");
    cfg.gap_extend          = cfg_file.Value("Scoring", "gap_extend");

    // Seeding parameters
    cfg.seed_shape_str      = (std::string) cfg_file.Value("Seed_params", "seed_shape");
    cfg.num_seeds_batch     = cfg_file.Value("Seed_params", "num_seeds_batch");
    cfg.chunk_size          = cfg_file.Value("Seed_params", "chunk_size");
    cfg.ignore_lower        = cfg_file.Value("Seed_params", "ignore_lower");
    cfg.use_transition      = cfg_file.Value("Seed_params", "use_transition");
    cfg.step                = cfg_file.Value("Seed_params", "step");

    // Filtering parameters
    cfg.xdrop               = cfg_file.Value("Filter_params", "xdrop");
    cfg.xdrop_threshold     = cfg_file.Value("Filter_params", "xdrop_threshold");

    // Extension parameters
    cfg.extension_threshold = cfg_file.Value("Extension_params", "extension_threshold");
    cfg.ydrop               = cfg_file.Value("Extension_params", "ydrop");
    cfg.do_gapped           = cfg_file.Value("Extension_params", "do_gapped");

    // Multi-threading
    cfg.num_threads         = 1;//tbb::task_scheduler_init::default_num_threads();

    //Output parameters
    cfg.output_format       = (std::string) cfg_file.Value("Output", "output_format");
    cfg.output_filename     = (std::string) cfg_file.Value("Output", "output_filename");

    gettimeofday(&end_time, NULL);

    useconds = end_time.tv_usec - start_time.tv_usec;
    seconds = end_time.tv_sec - start_time.tv_sec;
    mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    fprintf(stderr, "Time elapsed (loading configuration): %ld\n\n", mseconds);

    int nthreads = cfg.num_threads;
    tbb::task_scheduler_init init(nthreads);
    fprintf(stderr, "Using %d threads\n", cfg.num_threads);
    g_InitializeProcessor (0, 0);

    /////////// USER LOGIC ////////////////////
    g_DRAM = new DRAM;
    
    fprintf(stderr, "\nReading query file ...\n");
    gettimeofday(&start_time, NULL);
    gzFile f_rd = gzopen(cfg.query_filename.c_str(), "r");
    if (!f_rd) { fprintf(stderr, "cant open file: %s\n", cfg.query_filename.c_str()); exit(EXIT_FAILURE); }
        
    kseq_t *kseq_rd = kseq_init(f_rd);
    std::vector<seed_interval> interval_list;
    interval_list.clear();
    std::vector<uint32_t> chr_num_intervals;
    uint32_t q_chr_count = 0;
    uint32_t prev_num_intervals = 0;
    
    while (kseq_read(kseq_rd) >= 0) {
        size_t seq_len = kseq_rd->seq.l;
        std::string description = std::string(kseq_rd->name.s, kseq_rd->name.l);
        
        q_chr_coord.push_back(g_DRAM->bufferPosition);
        q_chr_id.push_back(description);
        q_chr_len.push_back(seq_len);

        if (g_DRAM->bufferPosition + seq_len > g_DRAM->size) {
            exit(EXIT_FAILURE); 
        }
        
        memcpy(g_DRAM->buffer + g_DRAM->bufferPosition, kseq_rd->seq.s, seq_len);
        g_DRAM->bufferPosition += seq_len;

        uint32_t curr_pos = 0;
        uint32_t end_pos = seq_len - 19;//sa->GetShapeSize();

        while (curr_pos < end_pos) {
            uint32_t start = curr_pos;
            uint32_t end = std::min(end_pos, start + cfg.num_seeds_batch);
            seed_interval inter;
            inter.start = start;
            inter.end = end;
            inter.num_invoked = 0;
            inter.num_intervals = 0;
            interval_list.push_back(inter);
            curr_pos += cfg.num_seeds_batch;
        }
        
        chr_num_intervals.push_back(interval_list.size()-prev_num_intervals);
        prev_num_intervals = interval_list.size();

        q_chr_count++;
    }
    g_DRAM->querySize = g_DRAM->bufferPosition;
    gzclose(f_rd);

    gettimeofday(&end_time, NULL);
    useconds = end_time.tv_usec - start_time.tv_usec;
    seconds = end_time.tv_sec - start_time.tv_sec;
    mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    fprintf(stderr, "Time elapsed (loading complete query from file): %ld msec \n\n", mseconds);

    f_rd = gzopen(cfg.reference_filename.c_str(), "r");
    if (!f_rd) { fprintf(stderr, "cant open file: %s\n", cfg.reference_filename.c_str()); exit(EXIT_FAILURE); }
        
    kseq_rd = kseq_init(f_rd);
    
    std::string q_chr;
    std::string r_chr;
    uint32_t q_len;
    uint32_t q_start;

    while (kseq_read(kseq_rd) >= 0) {

        short* chr_invoked = (short*) calloc(q_chr_count, sizeof(short));
        short* chr_sent = (short*) calloc(q_chr_count, sizeof(short));

        gettimeofday(&start_time, NULL);
        // reset bufferPosition to end of query
        g_DRAM->bufferPosition = g_DRAM->querySize;

        size_t seq_len = kseq_rd->seq.l;
        std::string description = std::string(kseq_rd->name.s, kseq_rd->name.l);

        fprintf(stderr, "Loading reference %s ...\n", description.c_str());

        if (g_DRAM->bufferPosition + seq_len > g_DRAM->size) {
            fprintf(stderr, "%ld exceeds DRAM size %ld \n", g_DRAM->bufferPosition + seq_len, g_DRAM->size);
            exit(EXIT_FAILURE); 
        }
        
        memcpy(g_DRAM->buffer + g_DRAM->bufferPosition, kseq_rd->seq.s, seq_len);

        sa = new SeedPosTable (g_DRAM->buffer, g_DRAM->bufferPosition, seq_len, cfg.seed_shape_str, cfg.step);

        g_SendRefWriteRequest (g_DRAM->bufferPosition, seq_len);
        g_DRAM->bufferPosition += seq_len;

        gettimeofday(&end_time, NULL);
        seconds = end_time.tv_sec - start_time.tv_sec;
        fprintf(stderr, "Time elapsed (loading reference and creating seed pos table): %ld sec \n\n", seconds);

        // start alignment
        tbb::flow::graph align_graph;

        printer_node printer(align_graph, tbb::flow::unlimited, segment_printer_body());

        tbb::flow::function_node<seeder_input, printer_input> seeder(align_graph, tbb::flow::unlimited, seeder_body());

        tbb::flow::make_edge(seeder, printer);

        tbb::flow::join_node<seeder_input> gatekeeper(align_graph);

        tbb::flow::make_edge(gatekeeper, seeder);

        tbb::flow::buffer_node<size_t> ticketer(align_graph);

        // Allocate tickets
        for (size_t t = 0ull; t < cfg.num_threads; t++)
            ticketer.try_put(t);

        tbb::flow::make_edge(tbb::flow::output_port<0>(printer), ticketer);

        tbb::flow::make_edge(ticketer, tbb::flow::input_port<1>(gatekeeper));

        uint32_t intervals_invoked = 0;
        uint32_t chr_intervals = 0;
        uint32_t num_intervals = 0;  
        uint32_t total_intervals = 0;  
        uint32_t prev_chr_intervals[2] = {0, 0};  
        uint32_t buffer = 0;
        bool new_chr = true;
        bond::blob chrom_seq;
        bond::blob chrom_rc_seq;
        char *rev_read_char = RevComp(chrom_seq);

        std::string send_q_chr;
        uint32_t send_q_len;
        uint32_t send_q_start;
        uint32_t send_buffer;
        uint32_t prev_buffer = 0;

        send_q_chr = q_chr_id[0];
        send_q_len = q_chr_len[0];
        send_q_start = q_chr_coord[0];
        send_buffer = 0;
        prev_buffer = 0;
        prev_chr_intervals[0] = chr_num_intervals[0]; 
        g_SendQueryWriteRequest (send_q_start, send_q_len, send_buffer);
        send_q_chr = q_chr_id[1];
        send_q_len = q_chr_len[1];
        send_q_start = q_chr_coord[1];
        send_buffer = 1;
        g_SendQueryWriteRequest (send_q_start, send_q_len, send_buffer);
        chr_sent[0] = 1;

        bool send_chr = false;
        bool invoke_chr = true; 
        uint32_t num_chr_sent = 2;
        uint32_t num_chr_invoked = 0;

        gettimeofday(&start_time, NULL);
        tbb::flow::source_node<seeder_payload> reader(align_graph,
            [&](seeder_payload &op) -> bool {
            while (true){
            if(num_chr_invoked < q_chr_count){  
                if(num_chr_invoked < num_chr_sent && invoke_chr){

                    q_chr = q_chr_id[num_chr_invoked];
                    q_len = q_chr_len[num_chr_invoked];
                    q_start = q_chr_coord[num_chr_invoked];
                    buffer = num_chr_invoked%2;
                    total_intervals += chr_intervals;
                    chr_intervals = chr_num_intervals[num_chr_invoked];

                    chrom_seq = bond::blob(g_DRAM->buffer + q_start, q_len);
                    rev_read_char = RevComp(chrom_seq);
                    chrom_rc_seq = bond::blob(rev_read_char, q_len);

                    intervals_invoked = 0;
                    invoke_chr = false;
                    fprintf(stderr, "Starting query %s %d ...\n", q_chr.c_str(), chr_intervals);
                }
                else if(num_chr_sent < q_chr_count && send_chr && prev_buffer == (num_chr_sent%2)){

                    send_q_chr = q_chr_id[num_chr_sent];
                    send_q_len = q_chr_len[num_chr_sent];
                    send_q_start = q_chr_coord[num_chr_sent];
                    prev_buffer = send_buffer;
                    prev_chr_intervals[prev_buffer] += chr_num_intervals[num_chr_sent-1];
                    send_buffer = num_chr_sent%2;

                    fprintf(stderr, "Sending query %s ...\n", send_q_chr.c_str());

                    g_SendQueryWriteRequest (send_q_start, send_q_len, send_buffer);
                    send_chr = false;
                    num_chr_sent++;
                }
                else{

                    fprintf(stderr, "else Starting query %s %d %d...\n", q_chr.c_str(), chr_intervals, intervals_invoked);
                    uint32_t curr_intervals_done;

                    if(prev_buffer == 0){
                        curr_intervals_done = seeder_body::num_seeded_regions0.load();
                    }
                    else{
                        curr_intervals_done = seeder_body::num_seeded_regions1.load();
                    }

                    if(chr_invoked > 0 && curr_intervals_done == prev_chr_intervals[prev_buffer]){
                        send_chr = true;
//                        fprintf(stderr, "Starting query %s %d ...\n", q_chr.c_str(), chr_intervals);
                    }

                    if (intervals_invoked < chr_intervals) {
                        seed_interval& inter = get<1>(op);
                        seed_interval curr_inter = interval_list[total_intervals + intervals_invoked++];
                        inter.start = curr_inter.start;
                        inter.end = curr_inter.end;
                        inter.num_invoked = intervals_invoked;
                        inter.num_intervals = chr_intervals;
                        inter.buffer = buffer;
                        reader_output& chrom = get<0>(op);
                        chrom.query_chr = q_chr;
                        chrom.ref_chr = description;
                        chrom.q_seq = chrom_seq;
                        chrom.q_rc_seq = chrom_rc_seq;
                        if(intervals_invoked == chr_intervals) {
                            invoke_chr = true;
                            num_chr_invoked++;
                        }
                        fprintf(stderr, "if Starting query %s %d %d...\n", q_chr.c_str(), chr_intervals, intervals_invoked);
                        return true;
                    }
                        return true;
                    fprintf(stderr, "end Starting query %s %d %d...\n", q_chr.c_str(), chr_intervals, intervals_invoked);
                }
            }
            else{
                gettimeofday(&start_time1, NULL);
                fprintf(stdout, "Completed\n");//Starting query ...%lu \n", start_time1.tv_sec);
                return false;
            }
            }
            printf("done %d\n", 
            }, true);

        tbb::flow::make_edge(reader, tbb::flow::input_port<0>(gatekeeper));
        
        align_graph.wait_for_all();

        gettimeofday(&end_time, NULL);
        seconds = end_time.tv_sec - start_time.tv_sec;
        fprintf(stderr, "Time elapsed (complete pipeline): %ld sec \n", seconds);
        delete[] rev_read_char;
    }

    gzclose(f_rd);
    
    fprintf(stderr, "#seeds: %lu \n", seeder_body::num_seeds.load());
    fprintf(stderr, "#seed hits: %lu \n", seeder_body::num_seed_hits.load());
    fprintf(stderr, "#anchors: %lu \n", seeder_body::num_hsps.load());

//    --------------------------------------------------------------------------
//     Shutdown and cleanup
//    -------------------------------------------------------------------------- 
    g_ShutdownProcessor();
}
