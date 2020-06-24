#include <iostream>
#include "parameters.h"

#define INVALID_KMER (1 << 31)
#define GRAIN_SIZE (1 << 18)

uint32_t NtChar2Int (char nt);
uint32_t TransitionNt (uint32_t nt);
void GenerateShapePos(std::string shape);
uint32_t GetKmerIndexAtPos(char* sequence, size_t pos, uint32_t seed_size);
int IsTransitionAtPos(int t);
