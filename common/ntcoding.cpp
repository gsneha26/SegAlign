#include <assert.h>
#include <string>
#include "ntcoding.h"
#include "parameters.h"

int shape_pos[32];
int shape_size;
int transition_pos[32];

inline uint32_t NtChar2Int (char nt) {
    switch(nt) {
        case 'A': return A_NT1;
        case 'C': return C_NT1;
        case 'G': return G_NT1;
        case 'T': return T_NT1;
        case 'N': return N_NT1;
        default: return N_NT1;
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
        if (nt[i] == N_NT1) {
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
