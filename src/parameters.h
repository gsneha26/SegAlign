#define TRANSITION_MASK 2
#define NUC 8 
#define NUC2 NUC*NUC
#define A_NT 0
#define C_NT 1
#define G_NT 2
#define T_NT 3
#define L_NT 4
#define N_NT 5
#define X_NT 6
#define E_NT 7

#define DEFAULT_LASTZ_INTERVAL 10000000
#define DEFAULT_WGA_CHUNK 250000

#define MAX_BLOCKS 1<<10
#define MAX_THREADS 1024 
#define BLOCK_SIZE 128 
#define NUM_WARPS 4
#define TILE_SIZE 32
#define BUFFER_DEPTH 1
#define SEQ_BLOCK_SIZE 500000000
