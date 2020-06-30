#include <assert.h>
#include "ntcoding.h"
#include "seed_filter.h"

class SeedPosTable {
    private:
        uint32_t index_table_size_;
        uint32_t ref_size_;
        int kmer_size_;
        int shape_size_;
        uint32_t *index_table_;
        uint32_t *pos_table_;
        uint32_t *tmp_pos_table_;

    public:
        SeedPosTable();
        SeedPosTable(char* ref_str, size_t start_addr, uint32_t ref_length, std::string shape, uint32_t step);
        ~SeedPosTable();

        int GetKmerSize();
        int GetShapeSize();
        int TouchKmerPos(std::string kmer); 
};
