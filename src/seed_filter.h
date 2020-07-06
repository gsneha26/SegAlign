#include "parameters.h"
#include "store.h"
#include <condition_variable>
#include <vector>

struct seedHit {
    uint32_t ref_start;
    uint32_t query_start;
};

struct hsp {
    uint32_t ref_start;
    uint32_t query_start;
    uint32_t len;
    int score;
};

typedef int(*InitializeProcessor_ptr)(int num_gpu);
typedef void(*InitializeSeeder_ptr)(bool transition, uint32_t WGA_CHUNK, uint32_t seed_size);
typedef void(*InitializeUngappedExtension_ptr)(int* sub_mat,int xdrop, int hspthresh, bool noentropy, int num_gpu);
typedef void(*SendSeedPosTable_ptr)(uint32_t* index_table, uint32_t index_table_size, uint32_t* pos_table, uint32_t ref_size, uint32_t max_pos_index);
typedef std::vector<hsp> (*SeedAndFilter_ptr)(std::vector<uint64_t> seed_offset_vector, bool rev, uint32_t buffer);
typedef uint32_t (*Filter_ptr)(char* r_seq, char* q_seq, uint32_t r_len, uint32_t q_len, uint32_t num_hits, seedHit* hits, hsp* hsp_out);
typedef void(*InclusivePrefixScan_ptr)(uint32_t* data, uint32_t len);
typedef void(*ShutdownUngappedExtension_ptr)();
typedef void(*ShutdownProcessor_ptr)();
typedef void(*clearQuery_ptr)(uint32_t buffer);
typedef void(*clearRef_ptr)();
typedef void(*SendRefWriteRequest_ptr)(size_t addr, size_t len);
typedef void(*SendQueryWriteRequest_ptr)(size_t addr, size_t len, uint32_t buffer);

extern InitializeProcessor_ptr g_InitializeProcessor;
extern InitializeSeeder_ptr g_InitializeSeeder;
extern InitializeUngappedExtension_ptr g_InitializeUngappedExtension;
extern SendSeedPosTable_ptr g_SendSeedPosTable;
extern SeedAndFilter_ptr g_SeedAndFilter;
extern Filter_ptr g_Filter;
extern SendRefWriteRequest_ptr g_SendRefWriteRequest;
extern SendQueryWriteRequest_ptr g_SendQueryWriteRequest;
extern InclusivePrefixScan_ptr g_InclusivePrefixScan;
extern ShutdownUngappedExtension_ptr g_ShutdownUngappedExtension;
extern ShutdownProcessor_ptr g_ShutdownProcessor;
extern clearQuery_ptr g_clearQuery;
extern clearRef_ptr g_clearRef;
