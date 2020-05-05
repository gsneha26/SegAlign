#include <string.h>
#include <boost/program_options.hpp> 
#include <boost/algorithm/string.hpp>
#include <tbb/task_scheduler_init.h>
#include <zlib.h>
#include <iostream>
#include "graph.h"
#include "kseq.h"
#include "store.h"

namespace po = boost::program_options;

////////////////////////////////////////////////////////////////////////////////
KSEQ_INIT2(, gzFile, gzread)

struct timeval start_time, end_time, start_time_complete, end_time_complete;
long useconds, seconds, mseconds;

Configuration cfg;
SeedPosTable *sa;

DRAM *ref_DRAM = nullptr;
DRAM *query_DRAM = nullptr;
DRAM *query_rc_DRAM = nullptr;

// query
std::vector<std::string> q_chr_id;
std::vector<uint32_t> q_buffer;
std::vector<uint32_t>  q_chr_len;
std::vector<std::string>  q_chr_index;
std::vector<size_t>  q_chr_coord;
std::vector<size_t>  rc_q_chr_coord;

// ref 
std::vector<std::string> r_chr_id;
std::vector<uint32_t>  r_chr_len;
std::vector<std::string>  r_chr_index;
std::vector<size_t>  r_chr_coord;

////////////////////////////////////////////////////////////////////////////////

char* RevComp(char* seq, uint32_t len) {

    char* rc = (char*) malloc(len*sizeof(char));

    for (size_t r = 0, i = len; i-- > 0;) {
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

//    gettimeofday(&start_time_complete, NULL);

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
        ("noentropy", po::bool_switch(&cfg.noentropy)->default_value(false), "don't use entropy factor for filtering stage")
        ("nounique", po::bool_switch(&cfg.nounique)->default_value(false), "don't remove repetitive elements")
        ("notrivial", po::bool_switch(&cfg.notrivial)->default_value(false), "Do not output a trivial self-alignment block if the target and query sequences are identical")
        ("output", po::value<std::string>(&cfg.output), "output filename")
        ("format", po::value<std::string>(&cfg.output_format)->default_value("maf-"), "format of output file")
        ("wga_chunk", po::value<uint32_t>(&cfg.wga_chunk_size)->default_value(DEFAULT_WGA_CHUNK), "chunk sizes for GPU calls for Xdrop - change only if you are a developer")
        ("lastz_interval", po::value<uint32_t>(&cfg.lastz_interval_size)->default_value(DEFAULT_LASTZ_INTERVAL), "LASTZ interval for ydrop - change only if you are a developer")
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

        fprintf(stderr, "Usage: run_wga_gpu target query data_folder [options]\n"); 
        std::cout << desc << std::endl;
        return 1;
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

    g_InitializeProcessor (cfg.sub_mat, cfg.transition, cfg.wga_chunk_size);

    ref_DRAM = new DRAM;
    query_DRAM = new DRAM;
    query_rc_DRAM = new DRAM;
    
    // Read query file
    if(cfg.debug){
	    gettimeofday(&start_time, NULL);
    }

    fprintf(stderr, "\nReading query file ...\n");

    gzFile f_rd = gzopen(cfg.query_filename.c_str(), "r");
    if (!f_rd) { 
        fprintf(stderr, "cant open file: %s\n", cfg.query_filename.c_str()); exit(EXIT_FAILURE); 
    }
        
    kseq_t *kseq_rd = kseq_init(f_rd);
    std::vector<seed_interval> interval_list;
    interval_list.clear();
    std::vector<uint32_t> chr_num_intervals;
    uint32_t total_q_chr = 0;
    uint32_t prev_num_intervals = 0;
    uint32_t total_query_intervals = 0;
    
    while (kseq_read(kseq_rd) >= 0) {
        size_t seq_len = kseq_rd->seq.l;
        std::string description = std::string(kseq_rd->name.s, kseq_rd->name.l);
        
        q_chr_coord.push_back(query_DRAM->bufferPosition);
        q_chr_id.push_back(description);
        q_buffer.push_back(0);
        q_chr_len.push_back(seq_len);
        q_chr_index.push_back(description);
        printf("%lu\n", query_DRAM->bufferPosition);

        if (query_DRAM->bufferPosition + seq_len > query_DRAM->size) {
            exit(EXIT_FAILURE); 
        }
        
        memcpy(query_DRAM->buffer + query_DRAM->bufferPosition, kseq_rd->seq.s, seq_len);
        query_DRAM->bufferPosition += seq_len;

        rc_q_chr_coord.push_back(query_rc_DRAM->bufferPosition);
        printf("%lu\n", query_rc_DRAM->bufferPosition);
        if (query_rc_DRAM->bufferPosition + seq_len > query_rc_DRAM->size) {
            exit(EXIT_FAILURE); 
        }
        char *query_rc = RevComp(kseq_rd->seq.s, seq_len);
        
        memcpy(query_rc_DRAM->buffer + query_rc_DRAM->bufferPosition, query_rc, seq_len);
        query_rc_DRAM->bufferPosition += seq_len;

        uint32_t curr_pos = 0;
        uint32_t end_pos = seq_len - cfg.seed_size;

        while (curr_pos < end_pos) {
            uint32_t start = curr_pos;
            uint32_t end = std::min(end_pos, start + cfg.lastz_interval_size);
            seed_interval inter;
            inter.start = start;
            inter.end = end;
            inter.num_invoked = 0;
            inter.num_intervals = 0;
            interval_list.push_back(inter);
            curr_pos += cfg.lastz_interval_size;
        }
        
        chr_num_intervals.push_back(interval_list.size()-prev_num_intervals);
        prev_num_intervals = interval_list.size();

        total_q_chr++;
    }

    total_query_intervals= interval_list.size();

    query_DRAM->seqSize = query_DRAM->bufferPosition;
    query_rc_DRAM->seqSize = query_rc_DRAM->bufferPosition;
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
    if (!f_rd) { 
        fprintf(stderr, "cant open file: %s\n", cfg.reference_filename.c_str()); exit(EXIT_FAILURE); 
    }
        
    kseq_rd = kseq_init(f_rd);
    uint32_t total_r_chr = 0;
    
    while (kseq_read(kseq_rd) >= 0) {
        size_t seq_len = kseq_rd->seq.l;
        std::string description = std::string(kseq_rd->name.s, kseq_rd->name.l);
        
        r_chr_coord.push_back(ref_DRAM->bufferPosition);
        r_chr_id.push_back(description);
        r_chr_len.push_back(seq_len);
        r_chr_index.push_back(description);

        printf("%lu\n", ref_DRAM->bufferPosition);
        if (ref_DRAM->bufferPosition + seq_len > ref_DRAM->size) {
            exit(EXIT_FAILURE); 
        }
        
        memcpy(ref_DRAM->buffer + ref_DRAM->bufferPosition, kseq_rd->seq.s, seq_len);
        ref_DRAM->bufferPosition += seq_len;

        total_r_chr++;
    }

    gzclose(f_rd);

    if(cfg.debug){
    	gettimeofday(&end_time, NULL);
    	useconds = end_time.tv_usec - start_time.tv_usec;
    	seconds = end_time.tv_sec - start_time.tv_sec;
    	mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    	fprintf(stderr, "Time elapsed (loading complete target from file): %ld msec \n\n", mseconds);
    }

    // start alignment
    fprintf(stderr, "\nStart alignment ...\n");
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

    std::string send_r_chr;
    uint32_t send_r_len;
    size_t send_r_start;
    bool send_ref_chr = true;
    uint32_t r_chr_sent = 0;

    uint32_t prev_chr_intervals[BUFFER_DEPTH];  
    uint32_t completed_intervals;  
    uint32_t chr_intervals_invoked;
    uint32_t chr_intervals_num;

    std::string send_q_chr;
    uint32_t send_q_len;
    size_t send_q_start;
    uint32_t send_buffer_id;
    bool send_query_chr = false;
    bool invoke_q_chr = false; 
    uint32_t q_chr_sent;
    std::string r_index;

    std::string r_chr;
    std::string q_chr;
    uint32_t q_len;
    size_t q_start;
    size_t rc_q_start;
    uint32_t q_buffer_id;
    uint32_t q_chr_invoked;
    std::string q_index;
    uint32_t q_num;
    uint32_t r_num;

    gettimeofday(&start_time, NULL);
    tbb::flow::source_node<seeder_payload> reader(align_graph,
        [&](seeder_payload &op) -> bool {

        while(true){
            if (send_ref_chr) {

                send_r_chr   = r_chr_id[r_chr_sent];
                send_r_len   = r_chr_len[r_chr_sent];
                send_r_start = r_chr_coord[r_chr_sent];
                r_index = r_chr_index[r_chr_sent];

                fprintf(stderr, "\nSending reference %s ...\n", send_r_chr.c_str());
                if(r_chr_sent > 0)
                    g_clearRef();
                g_SendRefWriteRequest (send_r_start, send_r_len);

                if(cfg.debug){
                    gettimeofday(&start_time_complete, NULL);
                }

                printf("%lu %lu\n", send_r_start, send_r_len);
                sa = new SeedPosTable (ref_DRAM->buffer, send_r_start, send_r_len, cfg.seed, cfg.step);
                if(cfg.debug){
                    gettimeofday(&end_time_complete, NULL);
                    useconds = end_time_complete.tv_usec - start_time_complete.tv_usec;
                    seconds = end_time_complete.tv_sec - start_time_complete.tv_sec;
                    mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
                    fprintf(stderr, "\nTime elapsed (seed position table create and copy to GPU) : %ld msec \n", mseconds);
                }

                for(int i = 0; i < BUFFER_DEPTH; i++){
                    seeder_body::num_seeded_regions[i] = 0;
                    prev_chr_intervals[i] = 0;
                }
                seeder_body::total_xdrop = 0;

                q_chr_invoked = 0;
                chr_intervals_invoked = 0;
                chr_intervals_num = 0;
                completed_intervals = 0;  
                q_buffer_id = 0;

                q_chr_sent = 0;
                send_query_chr = false;
                invoke_q_chr = true; 

                while(q_chr_sent < BUFFER_DEPTH && q_chr_sent < total_q_chr){

                    send_q_chr = q_chr_id[q_chr_sent];
                    send_q_len = q_chr_len[q_chr_sent];
                    send_q_start = q_chr_coord[q_chr_sent];
                    q_buffer[q_chr_sent] = q_chr_sent;
                    send_buffer_id = q_chr_sent;
                    fprintf(stderr, "\nSending query %s with buffer %d ...\n", send_q_chr.c_str(), send_buffer_id);
                    if(r_chr_sent > 0)
                        g_clearQuery(send_buffer_id);
                    g_SendQueryWriteRequest (send_q_start, send_q_len, send_buffer_id);
                    prev_chr_intervals[q_chr_sent] = chr_num_intervals[q_chr_sent];
                    q_chr_sent++;
                }

                r_chr_sent++;
                send_ref_chr = false;
            }
            else{
                if(q_chr_invoked > 0){
                    for(int i = 0; i < BUFFER_DEPTH; i++){
                        if(q_chr_sent < total_q_chr && seeder_body::num_seeded_regions[i] == prev_chr_intervals[i]){
                            send_q_chr = q_chr_id[q_chr_sent];
                            send_q_len = q_chr_len[q_chr_sent];
                            send_q_start = q_chr_coord[q_chr_sent];
                            q_buffer[q_chr_sent] = i;
                            prev_chr_intervals[i] += chr_num_intervals[q_chr_sent];

                            fprintf(stderr, "\nSending query %s with buffer %d ...\n", send_q_chr.c_str(), i);
                            g_clearQuery(i);
                            g_SendQueryWriteRequest (send_q_start, send_q_len, i);

                            q_chr_sent++;
                        }
                    }

                    if(r_chr_sent < total_r_chr && total_query_intervals == seeder_body::total_xdrop.load()){
                        send_ref_chr = true; 
                    }
                }
            }

            if(q_chr_invoked < q_chr_sent && invoke_q_chr) {

                q_chr = q_chr_id[q_chr_invoked];
                q_len = q_chr_len[q_chr_invoked];
                q_start = q_chr_coord[q_chr_invoked];
                rc_q_start = rc_q_chr_coord[q_chr_invoked];
                q_buffer_id = q_buffer[q_chr_invoked];
                q_index = q_chr_index[q_chr_invoked];
                completed_intervals += chr_intervals_num;
                chr_intervals_num = chr_num_intervals[q_chr_invoked];

                fprintf(stderr, "\nStarting query %s with buffer %d ...\n", q_chr.c_str(), q_buffer_id);

                chr_intervals_invoked = 0;
                invoke_q_chr = false;
            }

            if(q_chr_invoked < total_q_chr) {
                if (chr_intervals_invoked < chr_intervals_num) {
                    seed_interval& inter = get<1>(op);
                    seed_interval curr_inter = interval_list[completed_intervals + chr_intervals_invoked++];
                    inter.start = curr_inter.start;
                    inter.end = curr_inter.end;
                    inter.num_invoked = chr_intervals_invoked;
                    inter.num_intervals = chr_intervals_num;
                    inter.buffer = q_buffer_id;
                    reader_output& chrom = get<0>(op);
                    chrom.query_chr = q_chr;
                    chrom.ref_chr = send_r_chr;
                    chrom.q_start = q_start;
                    chrom.rc_q_start = rc_q_start;
                    chrom.q_len = q_len-cfg.seed_size;
                    chrom.q_index = q_index;
                    chrom.r_index = r_index;
                    if(chr_intervals_invoked == chr_intervals_num) {
                        q_chr_invoked++;
                        invoke_q_chr = true;
                    }
                    return true;
                }
            }
            else if(r_chr_sent == total_r_chr){
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

//    gettimeofday(&end_time_complete, NULL);
//    seconds = end_time_complete.tv_sec - start_time_complete.tv_sec;
//    fprintf(stderr, "\nTime elapsed (complete pipeline): %ld sec \n", seconds);
    return 0;
}
