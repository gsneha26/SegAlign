#include <condition_variable>
#include <vector>
#include "store.h"

/* Ungapped extension in this context uses the X-drop method for termination
 * A seedHit here refers to coordinates (ref, query) which is extended along the 
 * query forward strand (rev = 0) or along the reverse complement query strand
 * (rev = 1), It results in a segment (ref_start coordinate, query_start
 * coordinate, length of extension (left+right), score). Based on where or 
 * not the entropy should be factored (noentropy), the score of the segment is
 * adjusted. If the segment score crosses the threshold (hspthresh), the segment
 * qualifies as an HSP and the corresponding d_done is set to 1, otherwise 0.
 * Multiple seed hits could result in the same HSP. A prefix-sum of d_done helps
 * with gathering all the HSPs to the beginning. The resulting HSP vector is
 * sorted and the unique HSPs are generated. 
 */

// Each seed hit is 8B, segment is 16B
// With 32MB for the seed_hit array and both the HSPs array per 1GB GPU memory
// With higher GPU memory, the size just linearly increases

// seedHit vector is the input to the UngappedExtension function
struct seedHit {
    uint32_t ref_start;
    uint32_t query_start;
};

// format for the segment output which results after extending the seedHit
struct segment {
    uint32_t ref_start;
    uint32_t query_start;
    uint32_t len;
    int score;
};

// the function initializes the variables for ungapped xdrop extension
// num_gpu - number of GPUs that the process will be using, a d_hsp and d_done
// used during the extension is declared on each GPU
typedef void(*InitializeUngappedExtension_ptr)(int num_gpu, int* sub_mat, int input_xdrop, int input_hspthresh, bool input_noentropy);

// convert input sequence from alphabet to integers
typedef void(*CompressSeq_ptr)(char* input_seq, char* output_seq, uint32_t len);

// convert input sequence to its reverse complement and convert from alphabet to integers
typedef void(*CompressRevCompSeq_ptr)(char* input_seq, char* output_seq, uint32_t len);

//UngappedExtension function described above
//r_seq - compressed reference sequence (after CompressSeq() )
//q_seq - forward/reverse complement compressed query sequence
//r_len/q_len - length of r_seq/q_seq
//hits - seed hits that need to be extension
//num_hits - number of seed_hits
//hsp_out - HSPs that whose score crossed the threshold (size of the provided vector should be the same as num_hits) 
typedef uint32_t(*UngappedExtend_ptr)(char* r_seq, char* q_seq, uint32_t r_len, uint32_t q_len, uint32_t num_hits, seedHit* hits, segment* hsp_out);

// clear the vectors used in ungapped extension
typedef void(*ShutdownUngappedExtension_ptr)();

extern InitializeUngappedExtension_ptr g_InitializeUngappedExtension;
extern CompressSeq_ptr g_CompressSeq;
extern CompressRevCompSeq_ptr g_CompressRevCompSeq;
extern UngappedExtend_ptr g_UngappedExtend;
extern ShutdownUngappedExtension_ptr g_ShutdownUngappedExtension;

typedef int(*InitializeProcessor_ptr)(int num_gpu, bool transition, uint32_t WGA_CHUNK, uint32_t input_seed_size, int* sub_mat, int input_xdrop, int input_hspthresh, bool input_noentropy);
typedef void(*InclusivePrefixScan_ptr)(uint32_t* data, uint32_t len);
typedef void(*SendSeedPosTable_ptr)(uint32_t* index_table, uint32_t index_table_size, uint32_t* pos_table, uint32_t ref_size, uint32_t max_pos_index);
typedef void(*SendRefWriteRequest_ptr)(size_t addr, uint32_t len);
typedef void(*SendQueryWriteRequest_ptr)(size_t addr, uint32_t len, uint32_t buffer);
typedef std::vector<segment> (*SeedAndFilter_ptr)(std::vector<uint64_t> seed_offset_vector, bool rev, uint32_t buffer);
typedef void(*clearRef_ptr)();
typedef void(*clearQuery_ptr)(uint32_t buffer);
typedef void(*ShutdownProcessor_ptr)();

extern InitializeProcessor_ptr g_InitializeProcessor;
extern InclusivePrefixScan_ptr g_InclusivePrefixScan;
extern SendSeedPosTable_ptr g_SendSeedPosTable;
extern SendRefWriteRequest_ptr g_SendRefWriteRequest;
extern SendQueryWriteRequest_ptr g_SendQueryWriteRequest;
extern SeedAndFilter_ptr g_SeedAndFilter;
extern clearRef_ptr g_clearRef;
extern clearQuery_ptr g_clearQuery;
extern ShutdownProcessor_ptr g_ShutdownProcessor;
