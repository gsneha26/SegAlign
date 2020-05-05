#include "DRAM.h"
#include <vector>
#include <string>

extern DRAM *ref_DRAM;
extern DRAM *query_DRAM;
extern DRAM *query_rc_DRAM;
extern std::vector<std::string> q_chr_name;
extern std::vector<size_t>      q_chr_start;
extern std::vector<size_t>      q_chr_len;
extern std::vector<size_t>      q_chr_len_padded;
extern std::vector<std::string> r_chr_name;
extern std::vector<size_t>      r_chr_start;
extern std::vector<size_t>      r_chr_len;
extern std::vector<size_t>      r_chr_len_padded;
