#include <vector>
#include <string>
#include "DRAM.h"
#include "graph.h"

extern DRAM *ref_DRAM;
extern DRAM *query_DRAM;
extern DRAM *query_rc_DRAM;

extern std::vector<std::string> q_chr_name;
extern std::vector<uint32_t>    q_chr_file_name;
extern std::vector<size_t>      q_chr_start;
extern std::vector<uint32_t>    q_chr_len;
extern std::vector<std::string> rc_q_chr_name;
extern std::vector<uint32_t>    rc_q_chr_file_name;
extern std::vector<size_t>      rc_q_chr_start;
extern std::vector<uint32_t>    rc_q_chr_len;
extern std::vector<std::string> r_chr_name;
extern std::vector<uint32_t>    r_chr_file_name;
extern std::vector<size_t>      r_chr_start;
extern std::vector<uint32_t>    r_chr_len;
