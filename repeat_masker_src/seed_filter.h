#include <vector>
#include "graph.h"

typedef int(*InitializeProcessor_ptr)(int num_gpu, bool transition, uint32_t WGA_CHUNK, uint32_t input_seed_size, int* sub_mat, int input_xdrop, int input_hspthresh, bool input_noentropy);
typedef void(*SendRefWriteRequest_ptr)(char* seq, size_t addr, uint32_t len);
typedef std::vector<segment> (*SeedAndFilter_ptr)(std::vector<uint64_t> seed_offset_vector, bool rev, uint32_t ref_start, uint32_t ref_end);
typedef void(*clearRef_ptr)();
typedef void(*ShutdownProcessor_ptr)();

extern InitializeProcessor_ptr g_InitializeProcessor;
extern SendRefWriteRequest_ptr g_SendRefWriteRequest;
extern SeedAndFilter_ptr g_SeedAndFilter;
extern clearRef_ptr g_clearRef;
extern ShutdownProcessor_ptr g_ShutdownProcessor;
