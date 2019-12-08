#include "DRAM.h"                                                                                                                                                                        
#include <vector>
#include <string.h>
#include <stdio.h>
#include <mutex>

#define INF (1 << 28)
#define MAX_TILE_SIZE 512

#define MAX_NUM_TILES (1 << 20)
#define MASK_REF_REV 1
#define MASK_QUERY_REV 2
#define MASK_REF_COMP 3
#define MASK_QUERY_COMP 4

typedef int AlnOp;
enum AlnOperands { ZERO_OP, INSERT_OP, DELETE_OP, LONG_INSERT_OP, LONG_DELETE_OP};
enum states { Z, I, D, M , L_Z, L_I, L_D};

struct filter_tile {
    filter_tile (size_t ro, size_t qo, size_t rl, size_t ql, size_t qs)
        : ref_offset(ro),
          query_offset(qo),
          ref_length(rl),
          query_length(ql),
          query_tile_start(qs)
    {};
    size_t ref_offset;
    size_t query_offset;
    size_t ref_length;
    size_t query_length;
    size_t query_tile_start;
};

struct tile_output {
    tile_output (int id, int s, uint32_t ro, uint32_t qo)
      : batch_id(id),
        tile_score(s),
        max_ref_offset(ro),
        max_query_offset(qo)
    {};
    int batch_id;
    int tile_score;
    uint32_t max_ref_offset;
    uint32_t max_query_offset;

};

struct extend_tile {
    extend_tile (size_t ro, size_t qo, size_t rl, size_t ql) 
        : ref_offset(ro),
          query_offset(qo),
          ref_length(rl),
          query_length(ql)
    {};
    size_t ref_offset;
    size_t query_offset;
    size_t ref_length;
    size_t query_length;
};

struct extend_output {
    uint32_t max_ref_offset;
    uint32_t max_query_offset;
    std::vector<uint32_t> tb_pointers;

};

typedef size_t(*InitializeProcessor_ptr)(int t, int f);
typedef void(*SendRequest_ptr)(size_t ref_offset, size_t query_offset, size_t ref_length, size_t query_length, uint8_t align_fields);
typedef void(*SendSeedPosTable_ptr)(uint32_t* index_table, uint32_t index_table_size, uint64_t* pos_table, uint32_t ref_size);
typedef int (*SeedAndFilter_ptr)(std::vector<uint64_t> seed_offset_vector, int n);
typedef std::vector<tile_output> (*SendBatchRequest_ptr)(std::vector<filter_tile> tiles, uint8_t align_fields, int thresh);
typedef extend_output (*GACTXRequest_ptr)(extend_tile tile, uint8_t align_fields);
typedef void(*ShutdownProcessor_ptr)();
typedef void(*SendRefWriteRequest_ptr)(size_t addr, size_t len);
typedef void(*SendQueryWriteRequest_ptr)(size_t addr, size_t len);

extern DRAM *g_DRAM;
    
extern InitializeProcessor_ptr g_InitializeProcessor;
extern SendSeedPosTable_ptr g_SendSeedPosTable;
extern SeedAndFilter_ptr g_SeedAndFilter;
extern SendRequest_ptr g_SendRequest;
extern SendBatchRequest_ptr g_SendBatchRequest;
extern GACTXRequest_ptr g_GACTXRequest;
extern SendRefWriteRequest_ptr g_SendRefWriteRequest;
extern SendQueryWriteRequest_ptr g_SendQueryWriteRequest;
extern ShutdownProcessor_ptr g_ShutdownProcessor;      
