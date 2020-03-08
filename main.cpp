
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <boost/program_options.hpp> 
#include <boost/algorithm/string.hpp>
#include <tbb/task_scheduler_init.h>
#include <zlib.h>
#include <iostream>
#include "graph.h"
#include "kseq.h"
#include "DRAM.h"

namespace po = boost::program_options;

////////////////////////////////////////////////////////////////////////////////
KSEQ_INIT2(, gzFile, gzread)

struct timeval start_time, end_time, start_time_complete, end_time_complete;
long useconds, seconds, mseconds;

Configuration cfg;
SeedPosTable *sa;

// query
std::vector<std::string> q_chr_id;
std::vector<uint32_t>  q_chr_len;
std::vector<size_t>  q_chr_coord;

// ref 
std::vector<std::string> r_chr_id;
std::vector<uint32_t>  r_chr_len;
std::vector<size_t>  r_chr_coord;

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

    gettimeofday(&start_time_complete, NULL);

    po::options_description desc{"Options"};
    desc.add_options()
        ("scoring", po::value<std::string>(&cfg.scoring_file), "Scoring file in LASTZ format")
        ("ambiguous", po::value<std::string>(&cfg.ambiguous), "ambiguous nucleotides (n/iupac)")
        ("seed", po::value<std::string>(&cfg.seed_shape)->default_value("12of19"), "seed pattern")
        ("step", po::value<uint32_t>(&cfg.step)->default_value(1), "step length")
        ("xdrop", po::value<int>(&cfg.xdrop)->default_value(910), "x-drop threshold")
        ("ydrop", po::value<int>(&cfg.ydrop)->default_value(9430), "y-drop threshold")
        ("hspthresh", po::value<int>(&cfg.hspthresh)->default_value(3000), "threshold for high scoring pairs")
        ("gappedthresh", po::value<int>(&cfg.gappedthresh), "threshold for gapped alignments")
        ("notransition", po::bool_switch(&cfg.transition)->default_value(false), "allow (or don't) one transition in a seed hit")
        ("nogapped", po::bool_switch(&cfg.gapped)->default_value(false), "don't do gapped extension")
        ("format", po::value<std::string>(&cfg.output_format)->default_value("maf-"), "format of output file")
        ("debug", po::bool_switch(&cfg.debug)->default_value(false), "print debug messages")
        ("help", "Print help messages");

    po::options_description hidden;
    hidden.add_options()
        ("target", po::value<std::string>(&cfg.reference_filename)->required(), "target file")
        ("query", po::value<std::string>(&cfg.query_filename)->required(), "query file")
        ("data_folder", po::value<std::string>(&cfg.data_folder)->required(), "folder with 2bit files");

    po::options_description all_options;
    all_options.add(desc);
    all_options.add(hidden);

    po::positional_options_description p;
    p.add("target", 1);
    p.add("query", 1);
    p.add("data_folder", 1);

    po::variables_map vm;
    try{
        po::store(po::command_line_parser(argc, argv).options(all_options).positional(p).run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        if(!vm.count("help")){
            if(!vm.count("target") || !vm.count("query")){
		fprintf(stderr, "You must specify a target file and a query file\n"); 
            }
        }

        fprintf(stderr, "Usage: run_wga_gpu target query \"[options]\"\n"); 
        std::cout << desc << std::endl;
        return false;
    }

    cfg.transition = !cfg.transition;
    cfg.gapped = !cfg.gapped;
    if(cfg.seed_shape == "12of19"){
        cfg.seed = "TTT0T00TT00T0T0TTTT"; 
        cfg.seed_size = 19;
    }
    else if(cfg.seed_shape == "14of22"){
        cfg.seed = "TTT0T0TT00TT00T0T0TTTT";
        cfg.seed_size = 22;
    }
    else{
        int seed_len = cfg.seed_shape.size();
        cfg.seed = cfg.seed_shape;
        for(int i = 0; i< seed_len; i++){
            if(cfg.seed_shape[i] == '1')
                cfg.seed[i] = 'T';
            else
                cfg.seed[i] = '0';
        }
        cfg.seed_size = cfg.seed.size();
    }

    if(vm.count("gappedthresh") == 0)
        cfg.gappedthresh = cfg.hspthresh; 

    int ambiguous_reward = -100;
    int ambiguous_penalty = -100;
    int fill_score = -100;
    int bad_score = -1000;
    std::string ambiguous_field = "x";

    std::vector <std::string> fields;
    boost::split( fields, cfg.ambiguous, boost::is_any_of( "," ) );
    ambiguous_field = fields[0];
    if(fields.size() == 3){
        ambiguous_reward  = std::stoi(fields[1]);
        ambiguous_penalty = -1*std::stoi(fields[2]);
    }
    else if (cfg.ambiguous == "n" || cfg.ambiguous == "iupac"){
        ambiguous_reward  = 0;
        ambiguous_penalty = 0;
    }

    if(vm.count("scoring") == 0){

        //ACGT
        int tmp_sub_mat[L_NT][L_NT] = {{   91, -114,  -31, -123},
                                 { -114,  100, -125,  -31},
                                 {  -31, -125,  100, -114},
                                 { -123,  -31, -114,  91}};

        for(int i = 0; i < L_NT; i++){
            for(int j = 0; j < L_NT; j++){
                cfg.sub_mat[i*NUC+j] = tmp_sub_mat[i][j];
            }
        }

        //lower case characters
        for(int i = 0; i < L_NT; i++){
            cfg.sub_mat[i*NUC+L_NT] = bad_score;
            cfg.sub_mat[L_NT*NUC+i] = bad_score;
        }
        cfg.sub_mat[L_NT*NUC+L_NT] = bad_score;

        //N
        if(ambiguous_field == "n" || ambiguous_field == "iupac"){
            for(int i = 0; i < N_NT; i++){
                cfg.sub_mat[i*NUC+N_NT] = ambiguous_penalty;
                cfg.sub_mat[N_NT*NUC+i] = ambiguous_penalty;
            }
            cfg.sub_mat[N_NT*NUC+N_NT] = ambiguous_reward;
        }
        else{
            for(int i = 0; i < N_NT; i++){
                cfg.sub_mat[i*NUC+N_NT] = bad_score;
                cfg.sub_mat[N_NT*NUC+i] = bad_score;
            }
            cfg.sub_mat[N_NT*NUC+N_NT] = bad_score;
        }

        //other IUPAC
        if(ambiguous_field == "iupac"){
            for(int i = 0; i < X_NT; i++){
                cfg.sub_mat[i*NUC+X_NT] = ambiguous_penalty;
                cfg.sub_mat[X_NT*NUC+i] = ambiguous_penalty;
            }
            cfg.sub_mat[X_NT*NUC+X_NT] = ambiguous_reward;
        }
        else{
            for(int i = 0; i < L_NT; i++){
                cfg.sub_mat[i*NUC+X_NT] = fill_score;
                cfg.sub_mat[X_NT*NUC+i] = fill_score;
            }

            for(int i = L_NT; i < X_NT; i++){
                cfg.sub_mat[i*NUC+X_NT] = bad_score;
                cfg.sub_mat[X_NT*NUC+i] = bad_score;
            }
            cfg.sub_mat[X_NT*NUC+X_NT] = fill_score;
        }
    }

    cfg.num_threads = tbb::task_scheduler_init::default_num_threads();
    cfg.num_threads = (cfg.num_threads == 1) ? 2 : cfg.num_threads;
    tbb::task_scheduler_init init(cfg.num_threads);

    if(cfg.debug){
        fprintf(stderr, "Target %s\n", cfg.reference_filename.c_str());
        fprintf(stderr, "Query %s\n", cfg.query_filename.c_str());
        fprintf(stderr, "Seed %s\n", cfg.seed.c_str());
        fprintf(stderr, "Seed size %d\n", cfg.seed_size);
        fprintf(stderr, "Transition %d\n", cfg.transition);
        fprintf(stderr, "Gapped %d\n",cfg.gapped);
        fprintf(stderr, "xdrop %d\n", cfg.xdrop);
        fprintf(stderr, "ydrop %d\n", cfg.ydrop);
        fprintf(stderr, "HSP threshold %d\n", cfg.hspthresh);
        fprintf(stderr, "gapped threshold %d\n", cfg.gappedthresh);
        fprintf(stderr, "ambiguous %s\n", cfg.ambiguous.c_str());

        for(int i = 0; i < NUC; i++){
            for(int j = 0; j < NUC; j++){
                fprintf(stderr, "%d ", cfg.sub_mat[i*NUC+j]);
            }
            fprintf(stderr, "\n");
        }
    }

    fprintf(stderr, "Using %d threads\n", cfg.num_threads);

    g_InitializeProcessor (cfg.sub_mat);

    /////////// USER LOGIC ////////////////////
    g_DRAM = new DRAM;
    
    // Read query file
    if(cfg.debug){
	    gettimeofday(&start_time, NULL);
    }

    fprintf(stderr, "\nReading query file ...\n");

    gzFile f_rd = gzopen(cfg.query_filename.c_str(), "r");
    if (!f_rd) { fprintf(stderr, "cant open file: %s\n", cfg.query_filename.c_str()); exit(EXIT_FAILURE); }
        
    kseq_t *kseq_rd = kseq_init(f_rd);
    std::vector<seed_interval> interval_list;
    interval_list.clear();
    std::vector<uint32_t> chr_num_intervals;
    uint32_t q_chr_count = 0;
    uint32_t prev_num_intervals = 0;
    uint32_t total_query_intervals = 0;
    
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
        uint32_t end_pos = seq_len - cfg.seed_size;

        while (curr_pos < end_pos) {
            uint32_t start = curr_pos;
            uint32_t end = std::min(end_pos, start + LASTZ_INTERVAL);
            seed_interval inter;
            inter.start = start;
            inter.end = end;
            inter.num_invoked = 0;
            inter.num_intervals = 0;
            interval_list.push_back(inter);
            curr_pos += LASTZ_INTERVAL;
        }
        
        chr_num_intervals.push_back(interval_list.size()-prev_num_intervals);
        prev_num_intervals = interval_list.size();

        q_chr_count++;
    }

    total_query_intervals= interval_list.size();

    g_DRAM->querySize = g_DRAM->bufferPosition;
    gzclose(f_rd);

    if(cfg.debug){
    	gettimeofday(&end_time, NULL);
    	useconds = end_time.tv_usec - start_time.tv_usec;
    	seconds = end_time.tv_sec - start_time.tv_sec;
    	mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    	fprintf(stderr, "Time elapsed (loading complete query from file): %ld msec \n\n", mseconds);
    }

    // Read target file
    if(cfg.debug){
	    gettimeofday(&start_time, NULL);
    }

    fprintf(stderr, "\nReading target file ...\n");

    f_rd = gzopen(cfg.reference_filename.c_str(), "r");
    if (!f_rd) { fprintf(stderr, "cant open file: %s\n", cfg.reference_filename.c_str()); exit(EXIT_FAILURE); }
        
    kseq_rd = kseq_init(f_rd);
    uint32_t r_chr_count = 0;
    
    while (kseq_read(kseq_rd) >= 0) {
        size_t seq_len = kseq_rd->seq.l;
        std::string description = std::string(kseq_rd->name.s, kseq_rd->name.l);
        
        r_chr_coord.push_back(g_DRAM->bufferPosition);
        r_chr_id.push_back(description);
        r_chr_len.push_back(seq_len);

        if (g_DRAM->bufferPosition + seq_len > g_DRAM->size) {
            exit(EXIT_FAILURE); 
        }
        
        memcpy(g_DRAM->buffer + g_DRAM->bufferPosition, kseq_rd->seq.s, seq_len);
        g_DRAM->bufferPosition += seq_len;

        r_chr_count++;
    }

    gzclose(f_rd);

    if(cfg.debug){
    	gettimeofday(&end_time, NULL);
    	useconds = end_time.tv_usec - start_time.tv_usec;
    	seconds = end_time.tv_sec - start_time.tv_sec;
    	mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    	fprintf(stderr, "Time elapsed (loading complete target from file): %ld msec \n\n", mseconds);
    }

    fprintf(stderr, "\nStart alignment ...\n");
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

    std::string q_chr;
    std::string r_chr;
    uint32_t q_len;
    uint32_t q_start;
    bool invoke_ref_chr = true;

    uint32_t prev_chr_intervals[2] = {0, 0};  
    uint32_t intervals_invoked;
    uint32_t chr_intervals;
    uint32_t total_intervals;  
    uint32_t buffer;
    bond::blob chrom_seq;
    bond::blob chrom_rc_seq;
    char *rev_read_char = RevComp(chrom_seq);

    std::string send_q_chr;
    uint32_t send_q_len;
    uint32_t send_q_start;
    uint32_t send_buffer;
    std::string send_r_chr;
    uint32_t send_r_len;
    uint32_t send_r_start;
    uint32_t prev_buffer;
    uint32_t q_chr_sent;
    uint32_t r_chr_sent = 0;
    bool send_chr = false;
    bool invoke_chr = false; 
    uint32_t num_chr_invoked;
    std::string description;

    gettimeofday(&start_time, NULL);
    tbb::flow::source_node<seeder_payload> reader(align_graph,
        [&](seeder_payload &op) -> bool {

        while(true){
            if (invoke_ref_chr) {

                send_r_chr   = r_chr_id[r_chr_sent];
                send_r_len   = r_chr_len[r_chr_sent];
                send_r_start = r_chr_coord[r_chr_sent];

                fprintf(stderr, "\nSending reference %s ...\n", send_r_chr.c_str());
                if(r_chr_sent > 0)
                    g_clearRef();
                g_SendRefWriteRequest (send_r_start, send_r_len);
                sa = new SeedPosTable (g_DRAM->buffer, send_r_start, send_r_len, cfg.seed, cfg.step);

                seeder_body::num_seeded_regions0 = 0;
                seeder_body::num_seeded_regions1 = 0;
                intervals_invoked = 0;
                chr_intervals = 0;
                total_intervals = 0;  
                buffer = 0;
                prev_buffer = 0;
                q_chr_sent = 0;
                send_chr = false;
                invoke_chr = true; 
                num_chr_invoked = 0;

                while(q_chr_sent < 2 && q_chr_sent < q_chr_count){

                    send_q_chr = q_chr_id[q_chr_sent];
                    send_q_len = q_chr_len[q_chr_sent];
                    send_q_start = q_chr_coord[q_chr_sent];
                    send_buffer = q_chr_sent%2;
                    fprintf(stderr, "\nSending query %s ...\n", send_q_chr.c_str());
                    if(r_chr_sent > 0)
                        g_clearQuery(send_buffer);
                    g_SendQueryWriteRequest (send_q_start, send_q_len, send_buffer);
                    q_chr_sent++;
                }

                prev_chr_intervals[0] = chr_num_intervals[0]; 
                prev_chr_intervals[1] = 0;

                r_chr_sent++;
                invoke_ref_chr = false;

            }
            else if(send_chr && prev_buffer == (q_chr_sent%2)){

                send_q_chr = q_chr_id[q_chr_sent];
                send_q_len = q_chr_len[q_chr_sent];
                send_q_start = q_chr_coord[q_chr_sent];
                prev_buffer = send_buffer;
                prev_chr_intervals[prev_buffer] += chr_num_intervals[q_chr_sent-1];
                send_buffer = q_chr_sent%2;

                fprintf(stderr, "\nSending query %s ...\n", send_q_chr.c_str());
                g_clearQuery(send_buffer);
                g_SendQueryWriteRequest (send_q_start, send_q_len, send_buffer);

                send_chr = false;
                q_chr_sent++;
            }
            else{
                uint32_t curr_intervals_done;

                if(prev_buffer == 0){
                    curr_intervals_done = seeder_body::num_seeded_regions0.load();
                }
                else{
                    curr_intervals_done = seeder_body::num_seeded_regions1.load();
                }

                if(num_chr_invoked > 0 && curr_intervals_done == prev_chr_intervals[prev_buffer]){
                    if(q_chr_sent < q_chr_count)
                        send_chr = true;
                    else if(r_chr_sent < r_chr_count && total_query_intervals == (seeder_body::num_seeded_regions0.load()+seeder_body::num_seeded_regions1.load())){
                        invoke_ref_chr = true; 
                    }
                }
            }

            if(num_chr_invoked < q_chr_sent && invoke_chr) {

                q_chr = q_chr_id[num_chr_invoked];
                q_len = q_chr_len[num_chr_invoked];
                q_start = q_chr_coord[num_chr_invoked];
                buffer = num_chr_invoked%2;
                total_intervals += chr_intervals;
                chr_intervals = chr_num_intervals[num_chr_invoked];

                fprintf(stderr, "\nStarting query %s ...\n", q_chr.c_str());
                chrom_seq = bond::blob(g_DRAM->buffer + q_start, q_len);
                rev_read_char = RevComp(chrom_seq);
                chrom_rc_seq = bond::blob(rev_read_char, q_len);

                intervals_invoked = 0;
                invoke_chr = false;
            }

            if(num_chr_invoked < q_chr_count) {
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
                    chrom.ref_chr = send_r_chr;
                    chrom.q_seq = chrom_seq;
                    chrom.q_rc_seq = chrom_rc_seq;
                    if(intervals_invoked == chr_intervals) {
                        num_chr_invoked++;
                        if(num_chr_invoked < q_chr_count) {
                            invoke_chr = true;
                        }
                    }
                    return true;
                }
            }
            else if(r_chr_sent == r_chr_count){
                return false;
            }
            }

        }, true);

    tbb::flow::make_edge(reader, tbb::flow::input_port<0>(gatekeeper));
    
    align_graph.wait_for_all();

    if(cfg.debug){
        gettimeofday(&end_time, NULL);
        seconds = end_time.tv_sec - start_time.tv_sec;
        fprintf(stderr, "Time elapsed (complete pipeline): %ld sec \n", seconds);
    }

    delete[] rev_read_char;

    gzclose(f_rd);
    
    if(cfg.debug){
    	fprintf(stderr, "#seeds: %lu \n", seeder_body::num_seeds.load());
    	fprintf(stderr, "#seed hits: %lu \n", seeder_body::num_seed_hits.load());
    	fprintf(stderr, "#anchors: %lu \n", seeder_body::num_hsps.load());
    }

//    --------------------------------------------------------------------------
//     Shutdown and cleanup
//    -------------------------------------------------------------------------- 
    g_ShutdownProcessor();

    gettimeofday(&end_time_complete, NULL);
    seconds = end_time_complete.tv_sec - start_time_complete.tv_sec;
    fprintf(stderr, "\nTime elapsed (complete pipeline): %ld sec \n", seconds);
}
