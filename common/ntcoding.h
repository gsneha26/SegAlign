#define INVALID_KMER (1 << 31)
#define GRAIN_SIZE (1 << 18)

uint32_t NtChar2Int (char nt);
int GenerateShapePos(std::string shape);
int IsTransitionAtPos(int t);
uint32_t GetKmerIndexAtPos(char* sequence, size_t pos, uint32_t seed_size);
void RevComp(char* dst_buffer, char* src_buffer, size_t rc_start, size_t start, size_t len);
void GenerateSeedPosTable(char* ref_str, size_t start_addr, uint32_t ref_length, uint32_t step, int shape_size, int kmer_size);
