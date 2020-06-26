#pragma once
#define BOOST_LOCALE_NO_LIB
#define NOMINMAX
#include <tbb/flow_graph.h>
#include <tbb/reader_writer_lock.h>
#include <tbb/scalable_allocator.h>
#include <atomic>
#include <vector>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include "seed_pos_table.h"
#include "parameters.h"

using namespace tbb::flow;

struct Configuration {
    //FASTA files
    std::string reference_filename;
    std::string query_filename;
    std::string data_folder;
    std::string strand;

    //Scoring
    std::string scoring_file;
    std::string ambiguous;
    int sub_mat[NUC2];

    //Seed parameters
    std::string seed_shape;
    bool transition;
    uint32_t step;
    std::string seed;
    uint32_t seed_size;
    
    //Filter parameters
    int xdrop; 
    int hspthresh;
    bool noentropy;
    bool nounique;
    std::string hsplim;
    uint32_t hsplimit_num;

    //Extension parameters
    bool gapped;
    int ydrop;
    int gappedthresh;
    bool notrivial;

    // Output parameters
    std::string output_format;
    std::string output;
    bool markend;

    //System parameters
    uint32_t wga_chunk_size;
    uint32_t lastz_interval_size;
    int num_gpu;
    bool debug;

    uint32_t num_threads;
    uint32_t num_ref;
    uint32_t num_query;

};

extern Configuration cfg;
extern SeedPosTable *sa;

struct reader_output {
  size_t q_start;
  size_t r_start;
  size_t r_len;
  size_t q_len;
  uint32_t block_index;
  uint32_t r_block_index;
};

struct seed_interval {
    uint32_t start;
    uint32_t end;
    uint32_t num_invoked;
    uint32_t num_intervals;
    uint32_t buffer;
};

typedef std::vector<hsp> hsp_output; 

typedef tbb::flow::tuple <reader_output, seed_interval> seeder_payload;
typedef tbb::flow::tuple<int, hsp_output, hsp_output, uint32_t, size_t, size_t, size_t, size_t, size_t, uint32_t, uint32_t> printer_payload;
typedef tbb::flow::tuple <seeder_payload, size_t> seeder_input;
typedef tbb::flow::tuple<printer_payload, size_t> printer_input;

typedef tbb::flow::multifunction_node<printer_input, tbb::flow::tuple<size_t>> printer_node;

struct seeder_body{
	static std::atomic<uint64_t> num_seed_hits;
	static std::atomic<uint64_t> num_seeds;
	static std::atomic<uint64_t> num_hsps;
    	static std::atomic<uint32_t> total_xdrop;
    	static std::atomic<uint32_t> num_seeded_regions[BUFFER_DEPTH];
	printer_input operator()(seeder_input input);
};

struct segment_printer_body{
	void operator()(printer_input input, printer_node::output_ports_type & op);
};

