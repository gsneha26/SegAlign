#include <vector>
#include <cstdlib>
#include <stdint.h>
#include <math.h>
#include <string>
#include <assert.h>

struct seed_hit {
    uint32_t reference_offset;
    uint32_t query_offset;
};

struct Hits {
    Hits(uint64_t a, uint32_t b)
        : bin_offset(a),
        hit(b)
    {};

    uint64_t bin_offset;
    uint32_t hit;
};

static inline bool CompareHits (Hits h1, Hits h2) {
	return ((h1.bin_offset < h2.bin_offset));
}

class SeedPosTable {
    private:
        uint32_t index_table_size_;
        uint32_t ref_size_;
        int kmer_size_;
        int shape_size_;
        int bin_size_;
        uint32_t *index_table_;
        uint64_t *pos_table_;

    public:
        SeedPosTable();
        SeedPosTable(char* ref_str, uint32_t ref_length, std::string shape, int bin_size);
        ~SeedPosTable();

        int GetKmerSize();
        int GetShapeSize();
        std::vector<seed_hit> DSOFT(std::vector<uint64_t> seed_offset_vector, int threshold, uint32_t chunk_offset);
        int TouchKmerPos(std::string kmer); 
};

