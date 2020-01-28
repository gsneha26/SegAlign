#include "DRAM.h"                                                                                                                                                                        
#include <vector>
#include <string.h>
#include <stdio.h>
#include <mutex>

struct hsp {
    uint32_t ref_start;
    uint32_t query_start;
    uint32_t len;
    uint32_t score;
};

typedef size_t(*InitializeProcessor_ptr)(int t, int f);
typedef void(*SendSeedPosTable_ptr)(uint32_t* index_table, uint32_t index_table_size, uint64_t* pos_table, uint32_t ref_size);
typedef std::vector<hsp> (*SeedAndFilter_ptr)(std::vector<uint64_t> seed_offset_vector, bool rev);
typedef void(*ShutdownProcessor_ptr)();
typedef void(*SendRefWriteRequest_ptr)(size_t addr, size_t len);
typedef void(*SendQueryWriteRequest_ptr)(size_t addr, size_t len);

extern DRAM *g_DRAM;
    
extern InitializeProcessor_ptr g_InitializeProcessor;
extern SendSeedPosTable_ptr g_SendSeedPosTable;
extern SeedAndFilter_ptr g_SeedAndFilter;
extern SendRefWriteRequest_ptr g_SendRefWriteRequest;
extern SendQueryWriteRequest_ptr g_SendQueryWriteRequest;
extern ShutdownProcessor_ptr g_ShutdownProcessor;      
