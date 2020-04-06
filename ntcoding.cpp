#include "ntcoding.h"
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

inline uint32_t NtChar2IntCaseInsensitive (char nt) {
    switch(nt) {
        case 'a':
        case 'A': return A_NT;
        case 'c':
        case 'C': return C_NT;
        case 'g':
        case 'G': return G_NT;
        case 't':
        case 'T': return T_NT;
        case 'n':
        case 'N': return N_NT;
        default: return N_NT;
    }
}

uint32_t TransitionNt (uint32_t nt) {
    switch(nt) {
        case A_NT: return G_NT;
        case C_NT: return T_NT;
        case G_NT: return A_NT;
        case T_NT: return C_NT;
    }
}

void GenerateShapePos (std::string shape) {
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

int IsTransitionAtPos(int t) {
    return transition_pos[t];
}
