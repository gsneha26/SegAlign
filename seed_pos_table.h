#include <vector>
#include <cstdlib>
#include <stdint.h>
#include <math.h>
#include <string>
#include <assert.h>

class SeedPosTable {
    private:
        uint32_t index_table_size_;
        uint32_t ref_size_;
        int kmer_size_;
        int shape_size_;
        uint32_t *index_table_;
        uint64_t *pos_table_;
        uint32_t *tmp_pos_table_;

    public:
        SeedPosTable();
        SeedPosTable(char* ref_str, uint32_t ref_length, std::string shape);
        ~SeedPosTable();

        int GetKmerSize();
        int GetShapeSize();
        int TouchKmerPos(std::string kmer); 
};

