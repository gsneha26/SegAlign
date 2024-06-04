#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp> 
#include <iostream>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <tbb/reader_writer_lock.h>
#include <tbb/scalable_allocator.h>
#include <tbb/global_control.h>
#include "kseq.h"
#include "zlib.h"
#include "graph.h"
#include "ntcoding.h"
#include "seed_filter_interface.h"
#include "seed_filter.h"
#include "store.h"

namespace po = boost::program_options;

////////////////////////////////////////////////////////////////////////////////
KSEQ_INIT2(, gzFile, gzread)

struct timeval start_time, end_time, start_time_complete, end_time_complete;
long useconds, seconds, mseconds;

Configuration cfg;

DRAM *seq_DRAM = nullptr;
DRAM *seq_rc_DRAM = nullptr;

std::vector<std::string> chr_name;
std::vector<size_t>      chr_start;
std::vector<uint32_t>    chr_len;

std::vector<size_t>   block_start;
std::vector<uint32_t> block_len;
uint32_t total_chr = 0;

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv){

    po::options_description hidden;
    hidden.add_options()
        ("seq_file", po::value<std::string>(&cfg.seq_filename)->required(), "sequence file in FASTA format");

    po::options_description desc{"Sequence Options"};
    desc.add_options()
        ("strand", po::value<std::string>(&cfg.strand)->default_value("both"), "strand to search - plus/minus/both")
        ("neighbor_proportion", po::value<float>(&cfg.prop_neigh_interval)->default_value(0.2), "proportion of neighbouring intervals to align the query interval to");

    po::options_description scoring_desc{"Scoring Options"};
    scoring_desc.add_options()
        ("scoring", po::value<std::string>(&cfg.scoring_file), "Scoring file in LASTZ format")
        ("ambiguous", po::value<std::string>(&cfg.ambiguous), "ambiguous nucleotides - n/iupac");

    po::options_description seeding_desc{"Seeding Options"};
    seeding_desc.add_options()
        ("seed", po::value<std::string>(&cfg.seed_shape)->default_value("12of19"), "seed pattern-12of19(1110100110010101111)/14of22(1110101100110010101111)/an arbitrary pattern of 1s, 0s, and Ts ")
        ("step", po::value<uint32_t>(&cfg.step)->default_value(1), "Offset between the starting positions of successive target words considered for generating seed table")
        ("notransition", po::bool_switch(&cfg.seed.transition)->default_value(false), "don't allow one transition in a seed hit");

    po::options_description ungapped_desc{"Ungapped Extension Options"};
    ungapped_desc.add_options()
        ("xdrop", po::value<int>(&cfg.xdrop)->default_value(910), "x-drop value for gap-free extension")
        ("hspthresh", po::value<int>(&cfg.hspthresh)->default_value(3000), "segment score threshold for high scoring pairs")
        ("noentropy", po::bool_switch(&cfg.noentropy)->default_value(false), "don't adjust low score segment pair scores using entropy factor after filtering stage");

    po::options_description output_desc{"Output Options"};
    output_desc.add_options()
        ("M", po::value<uint32_t>(&cfg.M)->default_value(1), "report any position that is covered by at least this many alignments; the maximum allowed depth is 255")
        ("markend", po::bool_switch(&cfg.markend), "write a marker line just before completion");

    po::options_description system_desc{"System Options"};
    system_desc.add_options()
        ("wga_chunk_size", po::value<uint32_t>(&cfg.wga_chunk_size)->default_value(DEFAULT_WGA_CHUNK), "chunk sizes for GPU calls for Xdrop - change only if you are a developer")
        ("lastz_interval_size", po::value<uint32_t>(&cfg.lastz_interval_size)->default_value(DEFAULT_LASTZ_INTERVAL), "LASTZ interval for ydrop - change only if you are a developer")
        ("seq_block_size", po::value<uint32_t>(&cfg.seq_block_size)->default_value(DEFAULT_SEQ_BLOCK_SIZE), "LASTZ interval for ydrop - change only if you are a developer")
        ("num_gpu", po::value<int>(&cfg.num_gpu)->default_value(-1), "Specify number of GPUs to use - -1 if all the GPUs should be used")
        ("debug", po::bool_switch(&cfg.debug)->default_value(false), "print debug messages")
        ("version", "print version")
        ("help", "Print help messages");

    po::options_description all_options;
    all_options.add(desc);
    all_options.add(scoring_desc);
    all_options.add(seeding_desc);
    all_options.add(ungapped_desc);
    all_options.add(output_desc);
    all_options.add(system_desc);
    all_options.add(hidden);

    po::positional_options_description p;
    p.add("seq_file", 1);

    po::variables_map vm;
    try{
        po::store(po::command_line_parser(argc, argv).options(all_options).positional(p).run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        if(vm.count("version")){
            std::cerr << "SegAlign Version: " << VERSION << std::endl;
            return 0;
        }

        if(!vm.count("help")){
            if(!vm.count("seq_file")){
                fprintf(stderr, "You must specify a sequence file \n"); 
            }
        }

        fprintf(stderr, "\nUsage: segalign_repeat_masker seq_file [options]\n"); 
        fprintf(stderr, "Version: %s\n\n", VERSION);
        std::cerr << desc << std::endl;
        std::cerr << scoring_desc << std::endl;
        std::cerr << seeding_desc << std::endl;
        std::cerr << ungapped_desc << std::endl;
        std::cerr << output_desc << std::endl;
        std::cerr << system_desc << std::endl;
	
        if(vm.count("help"))
            return 0;
        else
            return 1;
    }

    cfg.seed.transition = !cfg.seed.transition;
    if(cfg.seed_shape == "12of19"){
        cfg.seed.shape = "TTT0T00TT00T0T0TTTT"; 
        cfg.seed.size = 19;
    }
    else if(cfg.seed_shape == "14of22"){
        cfg.seed.shape = "TTT0T0TT00TT00T0T0TTTT";
        cfg.seed.size = 22;
    }
    else{
        int seed_len = cfg.seed_shape.size();
        cfg.seed.shape = cfg.seed_shape;
        for(int i = 0; i< seed_len; i++){
            if(cfg.seed_shape[i] == '1')
                cfg.seed.shape[i] = 'T';
            else
                cfg.seed.shape[i] = '0';
        }
        cfg.seed.size = seed_len;
    }

    cfg.seed.kmer_size = GenerateShapePos(cfg.seed.shape);

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

        for(int i = 0; i < E_NT; i++){
            cfg.sub_mat[i*NUC+E_NT] = -10*cfg.xdrop;
            cfg.sub_mat[E_NT*NUC+i] = -10*cfg.xdrop;
        }
        cfg.sub_mat[E_NT*NUC+E_NT] = -10*cfg.xdrop;
    }

    tbb::global_control(tbb::global_control::max_allowed_parallelism, cfg.num_threads);

    if(cfg.debug){
        fprintf(stderr, "seq_file %s\n", cfg.seq_filename.c_str());
        fprintf(stderr, "ambiguous %s\n", cfg.ambiguous.c_str());
        fprintf(stderr, "Seed %s\n", cfg.seed.shape.c_str());
        fprintf(stderr, "Transition %d\n", cfg.seed.transition);
        fprintf(stderr, "xdrop %d\n", cfg.xdrop);
        fprintf(stderr, "HSP threshold %d\n", cfg.hspthresh);

        for(int i = 0; i < NUC; i++){
            for(int j = 0; j < NUC; j++){
                fprintf(stderr, "%d ", cfg.sub_mat[i*NUC+j]);
            }
            fprintf(stderr, "\n");
        }
    }

    fprintf(stderr, "Using %d threads\n", cfg.num_threads);

    cfg.num_gpu = g_InitializeInterface (cfg.num_gpu);
    g_InitializeProcessor (cfg.seed.transition, cfg.wga_chunk_size, cfg.seed.size, cfg.sub_mat, cfg.xdrop, cfg.hspthresh, cfg.noentropy);

    if(cfg.seq_block_size == DEFAULT_SEQ_BLOCK_SIZE){
        uint32_t offset = cfg.seq_block_size%cfg.lastz_interval_size;
        cfg.seq_block_size -= offset;
    }

    seq_DRAM = new DRAM;
    seq_rc_DRAM = new DRAM;
    
    gzFile f_rd;
    kseq_t *kseq_rd;

    // Read target file
    fprintf(stderr, "\nReading target file ...\n");

    if(cfg.debug){
	    gettimeofday(&start_time, NULL);
    }

    f_rd = gzopen(cfg.seq_filename.c_str(), "r");
    if (!f_rd) { 
        fprintf(stderr, "cant open file: %s\n", cfg.seq_filename.c_str()); 
        exit(7); 
    }
        
    kseq_rd = kseq_init(f_rd);

    while (kseq_read(kseq_rd) >= 0) {
        uint32_t seq_len = kseq_rd->seq.l;
        std::string seq_name = std::string(kseq_rd->name.s, kseq_rd->name.l);

        chr_name.push_back(seq_name);
        chr_start.push_back(seq_DRAM->bufferPosition);
        chr_len.push_back(seq_len);

        if (seq_DRAM->bufferPosition + seq_len > seq_DRAM->size) {
            fprintf(stderr, "Not enough memory in ref DRAM. Need to reduce input target sequence size to less than 6GB!\n"); 
            exit(9); 
        }
        
        memcpy(seq_DRAM->buffer + seq_DRAM->bufferPosition, kseq_rd->seq.s, seq_len);
        seq_DRAM->bufferPosition += seq_len;

        memset(seq_DRAM->buffer + seq_DRAM->bufferPosition, '&', 1);
        seq_DRAM->bufferPosition += 1;

        total_chr++;
    }

    seq_DRAM->bufferPosition -= 1;

    cfg.seq_len = seq_DRAM->bufferPosition;

    RevComp(seq_rc_DRAM->buffer, seq_DRAM->buffer, seq_rc_DRAM->bufferPosition, 0, cfg.seq_len);

    gzclose(f_rd);

    if(cfg.debug){
    	gettimeofday(&end_time, NULL);
    	useconds = end_time.tv_usec - start_time.tv_usec;
    	seconds = end_time.tv_sec - start_time.tv_sec;
    	mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    	fprintf(stderr, "Time elapsed (loading complete target from file): %ld msec \n\n", mseconds);
    }

    uint32_t total_query_intervals = ceil((float) cfg.seq_len/cfg.lastz_interval_size);
    cfg.num_neigh_interval = ceil((float) cfg.prop_neigh_interval*total_query_intervals);

    uint32_t left_intervals = ceil((float) (cfg.num_neigh_interval-1)/2); 
    uint32_t right_intervals = cfg.num_neigh_interval - 1 - left_intervals;
    uint32_t left_overlap = left_intervals * cfg.lastz_interval_size; 
    uint32_t right_overlap = right_intervals * cfg.lastz_interval_size; 
    uint32_t max_interval_seq_len = left_overlap + cfg.lastz_interval_size + right_overlap;

    if(cfg.debug)
        fprintf(stderr, "len: %lu lastz_interval: %u\ntotal_intervals: %u neigh_intervals: %u\nleft_intervals: %u 1 right_intervals: %u\nleft_overlap_size: %u\nright_overlap_size: %u\n\n", cfg.seq_len, cfg.lastz_interval_size, total_query_intervals, cfg.num_neigh_interval, left_intervals, right_intervals, left_overlap, right_overlap);

    size_t seq_block_start  = 0;
    uint32_t seq_block_len  = 0;
    uint32_t total_r_blocks = 0;

    std::vector<seed_interval> interval_list;
    interval_list.clear();
    std::vector<uint32_t> block_num_intervals;
    uint32_t prev_num_intervals = 0;
    uint32_t start_pos, end_pos;
    uint32_t interval_start, interval_end;
    bool left_overlap_limit = false;
    bool right_overlap_limit = false;

    for(size_t l = 0; l < cfg.seq_len; l+=cfg.seq_block_size){

        if(l < left_overlap) 
            seq_block_start = l;
        else
            seq_block_start = l-left_overlap;

        if(l+cfg.seq_block_size+right_overlap > cfg.seq_len)
            seq_block_len = static_cast<unsigned int>(cfg.seq_len - seq_block_start);
        else
            seq_block_len = static_cast<unsigned int>(l-seq_block_start+cfg.seq_block_size)+right_overlap;

        block_start.push_back(seq_block_start);
        block_len.push_back(seq_block_len);
        total_r_blocks += 1;

        start_pos = static_cast<unsigned int>(l-seq_block_start);
        if(seq_block_len < cfg.seq_block_size)
            end_pos = start_pos + seq_block_len - static_cast<unsigned int>(l-seq_block_start) - cfg.seed.size;
        else
            end_pos = start_pos + static_cast<unsigned int>(cfg.seq_block_size) - cfg.seed.size;

        if(cfg.debug)
            fprintf(stderr, "block l:%luM start:%luM len:%uM start_pos:%uM end_pos:%uM\n", l/1000000, seq_block_start/1000000, seq_block_len/1000000, start_pos/1000000, end_pos/1000000);

        while (start_pos < end_pos) {
            interval_start = start_pos;
            interval_end = std::min(end_pos, interval_start + cfg.lastz_interval_size);
            seed_interval inter;
            inter.start = interval_start;
            inter.end = interval_end;

            if(interval_start < left_overlap)
                left_overlap_limit = true;
            else
                left_overlap_limit = false;


            if((interval_end+right_overlap) > seq_block_len)
                right_overlap_limit = true;
            else
                right_overlap_limit = false;


            if(left_overlap_limit){
                inter.ref_start = 0;
                if(right_overlap_limit){
                    inter.ref_end = seq_block_len;
                }
                else{
                    if(max_interval_seq_len > seq_block_len)
                        inter.ref_end = seq_block_len;
                    else
                        inter.ref_end = max_interval_seq_len;
                }
            }
            else{
                if(right_overlap_limit){
                    inter.ref_end = seq_block_len;
                    if(seq_block_len < max_interval_seq_len)
                        inter.ref_start = 0;
                    else
                        inter.ref_start =  seq_block_len - max_interval_seq_len;
                }
                else{
                    inter.ref_start = interval_start-left_overlap;
                    inter.ref_end = interval_end+right_overlap;
                }
            }

            inter.num_invoked = 0;
            inter.num_intervals = 0;

            if(cfg.debug){
                float f1 = inter.start/1000000;
                float f2 = inter.end/1000000;
                fprintf(stderr, "%d %d | %u %u %u %u\n", left_overlap_limit, right_overlap_limit, inter.start, inter.end, inter.ref_start/1000000, inter.ref_end/1000000);
            }
            interval_list.push_back(inter);
            start_pos += cfg.lastz_interval_size;
        }

        block_num_intervals.push_back(interval_list.size()-prev_num_intervals);
        prev_num_intervals = interval_list.size();
    }

    total_query_intervals = interval_list.size();

//    for(int i = 0; i < 20; i++){
//        printf("%c %c\n", seq_DRAM->buffer[i], seq_rc_DRAM->buffer[cfg.seq_len-1-i]); 
//    }


    if(cfg.debug)
        for(int i = 0; i < total_r_blocks; i++)
            fprintf(stderr, "%u %lu %u\n", i, block_start[i]/1000000, block_len[i]/1000000);

    // start alignment
    fprintf(stderr, "\nStart alignment ...\n");
    tbb::flow::graph align_graph;

    printer_node printer(align_graph, tbb::flow::unlimited, interval_printer_body());

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

    size_t send_block_start;
    uint32_t send_block_len;

    uint32_t blocks_sent = 0;
    bool send_block = true;

    uint32_t prev_intervals_invoked = 0;
    uint32_t block_intervals_invoked = 0;
    uint32_t block_intervals_num = 0;

    seeder_body::total_xdrop = 0;

    gettimeofday(&start_time, NULL);
    tbb::flow::source_node<seeder_payload> reader(align_graph,
        [&](seeder_payload &op) -> bool {

        while(true){
            if (send_block) {

                send_block_start = block_start[blocks_sent];
                send_block_len   = block_len[blocks_sent];

                fprintf(stderr, "\nSending block %u ...\n", blocks_sent);

                if(blocks_sent > 0){
                    g_ClearRef();
                    g_ClearQuery();
                }

                g_SendRefWriteRequest (seq_DRAM->buffer, send_block_start, send_block_len);
                g_SendQueryWriteRequest ();

                if(cfg.debug){
                    gettimeofday(&start_time_complete, NULL);
                }

                GenerateSeedPosTable (seq_DRAM->buffer, send_block_start, send_block_len, cfg.step, cfg.seed.size, cfg.seed.kmer_size);

                if(cfg.debug){
                    gettimeofday(&end_time_complete, NULL);
                    useconds = end_time_complete.tv_usec - start_time_complete.tv_usec;
                    seconds = end_time_complete.tv_sec - start_time_complete.tv_sec;
                    mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
                    fprintf(stderr, "\nTime elapsed (seed position table create and copy to GPU) : %ld msec \n", mseconds);
                }

                seeder_body::num_seeded_regions = 0;

                send_block = false;

                prev_intervals_invoked += block_intervals_num;
                block_intervals_invoked = 0;
                block_intervals_num = block_num_intervals[blocks_sent];
            }
            else if(block_intervals_invoked == block_intervals_num && seeder_body::num_seeded_regions.load() == block_intervals_num){
                send_block = true;
            }

            if(blocks_sent < total_r_blocks) {
                if (block_intervals_invoked < block_intervals_num) {
                    seq_block& curr_block = std::get<0>(op);
                    curr_block.start = send_block_start;
                    curr_block.len   = send_block_len;
                    curr_block.index = blocks_sent; 

                    seed_interval& inter = std::get<1>(op);
                    seed_interval curr_inter = interval_list[prev_intervals_invoked + block_intervals_invoked++];
                    inter.start     = curr_inter.start;
                    inter.end       = curr_inter.end;
                    inter.ref_start = curr_inter.ref_start;
                    inter.ref_end   = curr_inter.ref_end;
                    inter.num_invoked   = block_intervals_invoked;
                    inter.num_intervals = block_intervals_num;

                    if(block_intervals_invoked == block_intervals_num) {
                        blocks_sent++;
                    }
                    return true;
                }
            }
            else {
                return false;
            }
            }

        }, true);

    tbb::flow::make_edge(reader, tbb::flow::input_port<0>(gatekeeper));
    
    align_graph.wait_for_all();

    g_ShutdownProcessor();

    if(cfg.debug){
        gettimeofday(&end_time, NULL);
        seconds = end_time.tv_sec - start_time.tv_sec;
        fprintf(stderr, "Time elapsed (complete pipeline): %ld sec \n\n", seconds);
    	fprintf(stderr, "#seeds: %lu \n", seeder_body::num_seeds.load());
    	fprintf(stderr, "#seed hits: %lu \n", seeder_body::num_seed_hits.load());
    	fprintf(stderr, "#HSPs: %lu \n", seeder_body::num_hsps.load());
    }

    return 0;
}
