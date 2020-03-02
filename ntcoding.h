#define A_NT 0
#define C_NT 1
#define G_NT 2
#define T_NT 3
#define L_NT 4
#define N_NT 5
#define X_NT 6

#include <stdint.h>
#include <iostream>
#include <assert.h>
#include <vector>

#define TRANSITION_MASK 2

uint32_t NtChar2Int (char nt);
uint32_t NtChar2IntCaseInsensitive (char nt);
uint32_t TransitionNt (uint32_t nt);
void GenerateShapePos(std::string shape);
uint32_t GetKmerIndexAtPos(char* sequence, uint32_t pos, uint32_t seed_size);
int IsTransitionAtPos(int t);
