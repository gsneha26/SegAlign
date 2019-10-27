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

#define do_traceback (1 << 5)
#define reverse_ref (1 << 4)
#define complement_ref (1 << 3)
#define reverse_query (1 << 2)
#define complement_query (1 << 1)
#define start_end 1

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

	// D-SOFT parameters
    std::string seed_shape_str;
	int bin_size;
	int dsoft_threshold;
	int num_seeds_batch;
    int chunk_size;
	int seed_occurence_multiple;
	int max_candidates;
	int num_nz_bins;
    bool ignore_lower;
    bool use_transition;
    int hash_size;
    
	// GACT scoring
	int gact_sub_mat[11];
	int gap_open;
	int gap_extend;

    // Banded GACT filter
    int xdrop; 
    int xdrop_limit;
    int xdrop_score_threshold;

    // GACT-X
    int tile_size;
    int tile_overlap;
    int extension_threshold;
    int ydrop;

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

typedef tbb::flow::tuple <reader_output, seed_interval> seeder_payload;
typedef tbb::flow::tuple <seeder_payload, size_t> seeder_input;

struct seeder_output {
    std::vector<seed_hit> fwHits;
    std::vector<seed_hit> rcHits;
};

typedef tbb::flow::tuple<reader_output, seeder_output> filter_payload;
typedef tbb::flow::tuple<filter_payload, size_t> filter_input;

struct anchor {
    anchor(uint32_t ro, uint32_t qo, int score) 
        : reference_offset(ro),
        query_offset(qo),
        score(score)
    {};
    uint32_t reference_offset;
    uint32_t query_offset;
    int score;
};

static inline bool CompareAnchors (anchor a1, anchor a2) {
	return ((a1.score > a2.score) || ((a1.score == a2.score) && (a1.reference_offset < a2.reference_offset)) || ((a1.score == a2.score) && (a1.reference_offset == a2.reference_offset) && (a1.query_offset == a2.query_offset)));
}

typedef std::vector<anchor> filter_output;

typedef tbb::flow::tuple<reader_output, filter_output, filter_output> extender_payload;
typedef tbb::flow::tuple<extender_payload, size_t> extender_input;

struct Alignment {
	int chr_id;
	uint32_t curr_reference_offset;
	uint32_t curr_query_offset;
	uint32_t reference_start_offset;
	uint32_t query_start_offset;
	uint32_t reference_end_offset;
	uint32_t query_end_offset;
	uint32_t reference_start_addr;
	uint32_t query_start_addr;
	uint32_t reference_length;
	uint32_t query_length;
	std::string aligned_reference_str;
	std::string aligned_query_str;
    int score;
	char strand;

    int num_left_tiles;
    int num_right_tiles;
};

typedef std::vector<Alignment> extender_output;

typedef tbb::flow::tuple<reader_output, extender_output> printer_payload;
typedef tbb::flow::tuple<printer_payload, size_t> printer_input;

struct seeder_body
{
	static std::atomic<uint64_t> num_seed_hits;
	static std::atomic<uint64_t> num_seeds;
	filter_input operator()(seeder_input input);
};

struct filter_body 
{
	static std::atomic<uint64_t> num_filter_tiles;
	static std::atomic<uint64_t> num_anchors;
	extender_input operator()(filter_input input);
};

typedef tbb::flow::multifunction_node<extender_input, tbb::flow::tuple<printer_input, size_t>> extender_node;

struct extender_body
{
	static std::atomic<uint64_t> num_extend_tiles;
	void operator()(extender_input input, extender_node::output_ports_type & op);
	Alignment makeAlignment(reader_output read, anchor anc, char strand);
};

struct maf_printer_body
{
	size_t operator()(printer_input input);
};
