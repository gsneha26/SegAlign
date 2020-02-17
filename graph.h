#pragma once
#define BOOST_LOCALE_NO_LIB
#define NOMINMAX
#include <algorithm>

#include <tbb/flow_graph.h>
#include <tbb/reader_writer_lock.h>
#include <tbb/scalable_allocator.h>
#include <mutex>
#include <atomic>
#include <vector>
#include <stdio.h>
#include <string>
#include <cstddef>

#include <bond/core/blob.h>
#include <bond/core/reflection.h>
#include "GPU.h"
#include "ConfigFile.h"
#include "seed_pos_table.h"

using namespace tbb::flow;

//reference
extern std::vector<std::string> r_chr_id;
extern std::vector<uint32_t>  r_chr_len;
extern std::vector<uint32_t>  r_chr_len_unpadded;
extern std::vector<uint32_t>  r_chr_coord;

extern FILE *mafFile;

struct Configuration {
    //FASTA files
    std::string reference_name;
    std::string query_name;
    std::string reference_filename;
    std::string query_filename;
    std::string data_folder;

    //Scoring
    int gact_sub_mat[11];
    int gap_open;
    int gap_extend;

    //Seed parameters
    std::string seed_shape_str;
    int num_seeds_batch;
    int chunk_size;
    bool ignore_lower;
    bool use_transition;
    
    //Filter parameters
    int xdrop; 
    int xdrop_threshold;

    //Extension parameters
    int extension_threshold;
    int ydrop;
    std::string lastz_path;
    bool do_gapped;

    //Multi-threading
    int num_threads;

    // Output
    std::string output_filename;
};

extern Configuration cfg;
extern SeedPosTable *sa;

struct reader_output {
	std::string description;
	bond::blob seq;
	bond::blob rc_seq;
};

struct seed_interval {
    uint32_t start;
    uint32_t end;
    uint32_t num_invoked;
    uint32_t num_intervals;
};

typedef std::vector<hsp> hsp_output; 

typedef tbb::flow::tuple <reader_output, seed_interval> seeder_payload;
typedef tbb::flow::tuple<int, hsp_output, hsp_output, std::string> printer_payload;
typedef tbb::flow::tuple <seeder_payload, size_t> seeder_input;
typedef tbb::flow::tuple<printer_payload, size_t> printer_input;

typedef tbb::flow::multifunction_node<printer_input, tbb::flow::tuple<size_t>> printer_node;

struct seeder_body{
	static std::atomic<uint64_t> num_seed_hits;
	static std::atomic<uint64_t> num_seeds;
	static std::atomic<uint64_t> num_hsps;
	printer_input operator()(seeder_input input);
};

struct segment_printer_body{
	void operator()(printer_input input, printer_node::output_ports_type & op);
};

