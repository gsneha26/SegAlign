#include <vector>
#include "graph.h"

typedef void(*InitializeProcessor_ptr)(bool transition, uint32_t WGA_CHUNK, uint32_t input_seed_size, int* sub_mat, int input_xdrop, int input_hspthresh, bool input_noentropy);
typedef void(*SendQueryWriteRequest_ptr)();
typedef std::vector<segmentPair> (*SeedAndFilter_ptr)(std::vector<uint64_t> seed_offset_vector, bool rev, uint32_t ref_start, uint32_t ref_end);
typedef void(*ClearQuery_ptr)();
typedef void(*ShutdownProcessor_ptr)();

extern InitializeProcessor_ptr g_InitializeProcessor;
extern SendQueryWriteRequest_ptr g_SendQueryWriteRequest;
extern SeedAndFilter_ptr g_SeedAndFilter;
extern ClearQuery_ptr g_ClearQuery;
extern ShutdownProcessor_ptr g_ShutdownProcessor;
