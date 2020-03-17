#include <iostream>
#include "parameters.h"

uint32_t NtChar2Int (char nt);
uint32_t NtChar2IntCaseInsensitive (char nt);
uint32_t TransitionNt (uint32_t nt);
void GenerateShapePos(std::string shape);
uint32_t GetKmerIndexAtPos(char* sequence, size_t pos, uint32_t seed_size);
int IsTransitionAtPos(int t);
