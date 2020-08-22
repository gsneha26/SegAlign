#include <assert.h>
#include <string>
#include "ntcoding.h"
#include "parameters.h"
#include "seed_filter.h"
#include "tbb/parallel_sort.h"

int shape_pos[32];
int shape_size;
int transition_pos[32];

inline uint32_t NtChar2Int (char nt) {
    switch(nt) {
        case 'A': return A_NT;
        case 'C': return C_NT;
        case 'G': return G_NT;
        case 'T': return T_NT;
        case 'N': return N_NT;
        default: return N_NT;
    }
}

int GenerateShapePos (std::string shape) {
    shape_size = 0;
    int j = 0;
    for (int i = 0; i < shape.length(); i++) {
        if ((shape[i] == '1') || (shape[i] == 'T')) {
            shape_pos[shape_size++] = i;
            if (shape[i] == 'T') {
                transition_pos[j] = 1;
            }
            else {
                transition_pos[j] = 0;
            }
            j++;
        }
    }
    return shape_size;
}

int IsTransitionAtPos(int t) {
    return transition_pos[t];
}

uint32_t GetKmerIndexAtPos (char* sequence, size_t pos, uint32_t seed_size) {

    uint32_t nt[64];

    for(int i = 0; i < seed_size; i++){
        nt[i] = NtChar2Int(sequence[pos+i]);
        if (nt[i] == N_NT) {
            return INVALID_KMER;
        }
    }

    uint32_t kmer = 0;

    for (int i = 0; i < shape_size; i++) {
        kmer = (kmer << 2) + nt[shape_pos[i]];
    }

    return kmer;
}

void RevComp(char* dst_buffer, char* src_buffer, size_t rc_start, size_t start, size_t len){

    size_t r = rc_start;
    for (size_t i = start+len; i> start; i--) {
        
        switch (src_buffer[i-1]) {
            case 'a': dst_buffer[r++] = 't';
                      break;

            case 'A': dst_buffer[r++] = 'T';
                      break;

            case 'c': dst_buffer[r++] = 'g';
                      break;

            case 'C': dst_buffer[r++] = 'G';
                      break;

            case 'g': dst_buffer[r++] = 'c';
                      break;

            case 'G': dst_buffer[r++] = 'C';
                      break;

            case 't': dst_buffer[r++] = 'a';
                      break;

            case 'T': dst_buffer[r++] = 'A';
                      break;

            case 'n': dst_buffer[r++] = 'n';
                      break;

            case 'N': dst_buffer[r++] = 'N';
                      break;

            case '&': dst_buffer[r++] = '&';
                      break;

            default: printf("Bad Nt char! '%c' %lu\n", src_buffer[i], i);
        }
    }
}

void GenerateSeedPosTable(char* ref_str, size_t start_addr, uint32_t ref_length, uint32_t step, int shape_size, int kmer_size) {

    assert(kmer_size <= 15);
    assert(kmer_size > 3); 

    uint32_t *index_table;
    uint32_t *pos_table;
    uint32_t index_table_size;

    uint32_t offset = (shape_size+1)%step;
    uint32_t start_offset = step - offset;
    
    index_table_size = ((uint32_t)1 << 2*kmer_size) + 1;
    index_table = (uint32_t*) calloc(index_table_size, sizeof(uint32_t));

    uint32_t num_steps = (ref_length - shape_size + offset) / step;

    uint32_t* tmp_index_arr = (uint32_t*) malloc(num_steps * sizeof(uint32_t));
    uint32_t* tmp_off_arr = (uint32_t*) malloc(num_steps * sizeof(uint32_t));

    tbb::parallel_for( tbb::blocked_range<uint32_t>(0, num_steps, GRAIN_SIZE),
            [&](tbb::blocked_range<uint32_t> r){

            for (uint32_t i=r.begin(); i<r.end(); ++i){
                uint32_t index = GetKmerIndexAtPos(ref_str, start_addr+start_offset+(i*step), shape_size); 
                tmp_index_arr[i] = index;

                // valid index
                if (index != (uint32_t) INVALID_KMER) {
                    tmp_off_arr[i] = __sync_fetch_and_add( &index_table[index+1], 1);
                }
            }
        });

    g_InclusivePrefixScan(index_table, index_table_size);

    uint32_t num_index = index_table[index_table_size-1];
    
    pos_table = (uint32_t*) malloc(num_index * sizeof(uint32_t));

    tbb::parallel_for( tbb::blocked_range<uint32_t>(0, num_steps, GRAIN_SIZE),
            [&](tbb::blocked_range<uint32_t> r){

            for (uint32_t i=r.begin(); i<r.end(); ++i){
                uint32_t index = tmp_index_arr[i];

                // valid index
                if (index != (uint32_t) INVALID_KMER) {
                    uint32_t curr_idx = index_table[index] + tmp_off_arr[i]; 
                    pos_table[curr_idx] = start_offset+(i*step);
                }       
            }
        });

    g_SendSeedPosTable(index_table+1, index_table_size-1, pos_table, num_index);

    free(tmp_index_arr);
    free(tmp_off_arr);
    free(index_table);
    free(pos_table);
}
