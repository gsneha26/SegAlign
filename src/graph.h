#pragma once
#define BOOST_LOCALE_NO_LIB
#define NOMINMAX
#include <atomic>
#include <string>
#include <tbb/flow_graph.h>
#include <vector>
#include "parameters.h"

#define DEFAULT_SEQ_BLOCK_SIZE 500000000
#define DEFAULT_LASTZ_INTERVAL 10000000
#define DEFAULT_WGA_CHUNK 250000

#define BUFFER_DEPTH 2

typedef struct Seed_config {
    std::string shape;
    int size;
    int kmer_size;
    bool transition;
} Seed_config;

typedef struct segmentPair {
    uint32_t ref_start;
    uint32_t query_start;
    uint32_t len;
    int score;
} segmentPair;

typedef struct Configuration {
    //Input files/folder name
    std::string reference_filename;
    std::string query_filename;
    std::string data_folder;

    // Sequence parameters
    std::string strand;

    // Scoring
    std::string scoring_file;
    std::string ambiguous;
    int sub_mat[NUC2];

    // Seed parameters
    Seed_config seed;
    std::string seed_shape;
    uint32_t step;
    
    // Filter parameters
    int xdrop; 
    int hspthresh;
    bool noentropy;

    // Extension parameters
    bool gapped;
    int ydrop;
    int gappedthresh;
    bool notrivial;

    // Output parameters
    std::string output_format;
    std::string output;
    bool markend;

    // System parameters
    uint32_t wga_chunk_size;
    uint32_t lastz_interval_size;
    uint32_t seq_block_size;
    int num_gpu;
    int num_threads;
    bool debug;
} Configuration;

extern Configuration cfg;

typedef struct seq_block {
  int r_index;
  int q_index;
  size_t r_start;
  size_t q_start;
  uint32_t r_len;
  uint32_t q_len;
} seq_block;

typedef struct seed_interval {
    uint32_t start;
    uint32_t end;
    uint32_t num_invoked;
    uint32_t num_intervals;
    uint32_t buffer;
} seed_interval;

typedef std::vector<segmentPair> hsp_output; 
typedef std::tuple<seq_block, seed_interval> seeder_payload;
typedef std::tuple<seeder_payload, hsp_output, hsp_output> printer_payload;
typedef std::tuple<seeder_payload, size_t> seeder_input;
typedef std::tuple<printer_payload, size_t> printer_input;
typedef tbb::flow::multifunction_node<printer_input, std::tuple<size_t>> printer_node;

typedef struct seeder_body {
	static std::atomic<uint64_t> num_seed_hits;
	static std::atomic<uint64_t> num_seeds;
	static std::atomic<uint64_t> num_hsps;
    static std::atomic<uint32_t> total_xdrop;
    static std::atomic<uint32_t> num_seeded_regions[BUFFER_DEPTH];
	printer_input operator()(seeder_input input);
} seeder_body;

typedef struct segment_printer_body {
	void operator()(printer_input input, printer_node::output_ports_type & op);
} segment_printer_body;

