#include "cuda_utils.h"
#include "parameters.h"
#include "seed_filter_interface.h"
#include "store_gpu.h"

#include <claraparabricks/genomeworks/cudaextender/extender.hpp>
#include <claraparabricks/genomeworks/cudaextender/utils.hpp>
#include <claraparabricks/genomeworks/utils/pinned_host_vector.hpp>
#include <claraparabricks/genomeworks/utils/cudautils.hpp>
using namespace claraparabricks::genomeworks;
using namespace cudaextender;
using namespace cudautils;

// Control Variables
std::mutex mu;
std::condition_variable cv;
std::vector<int> available_gpus;

int NUM_DEVICES;
int8_t** d_ref_seq;
uint32_t ref_len;

uint32_t** d_index_table;
uint32_t** d_pos_table;

__global__
void compress_string (uint32_t len, char* src_seq, char* dst_seq){ 
    int thread_id = threadIdx.x;
    int block_dim = blockDim.x;
    int grid_dim = gridDim.x;
    int block_id = blockIdx.x;

    int stride = block_dim * grid_dim;
    uint32_t start = block_dim * block_id + thread_id;

    for (uint32_t i = start; i < len; i += stride) {
        char ch = src_seq[i];
        char dst = X_NT1;
        if (ch == 'A')
            dst = A_NT1;
        else if (ch == 'C')
            dst = C_NT1;
        else if (ch == 'G')
            dst = G_NT1;
        else if (ch == 'T')
            dst = T_NT1;
        else if ((ch == 'a') || (ch == 'c') || (ch == 'g') || (ch == 't'))
            dst = L_NT1;
        else if ((ch == 'n') || (ch == 'N'))
            dst = N_NT1;
        else if (ch == '&')
            dst = E_NT1;
        dst_seq[i] = dst;
    }
}

int InitializeInterface (int num_gpu){

    int nDevices;

    cudaError_t err = cudaGetDeviceCount(&nDevices);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: No GPU device found!\n");
        exit(1);
    }

    if(num_gpu == -1){
        NUM_DEVICES = nDevices; 
    }
    else{
        if(num_gpu <= nDevices){
            NUM_DEVICES = num_gpu;
        }
        else{
            fprintf(stderr, "Requested GPUs greater than available GPUs\n");
            exit(10);
        }
    }

    fprintf(stderr, "Using %d GPU(s)\n", NUM_DEVICES);

    d_ref_seq = (int8_t**) malloc(NUM_DEVICES*sizeof(int8_t*));
    
    d_index_table = (uint32_t**) malloc(NUM_DEVICES*sizeof(uint32_t*));
    d_pos_table = (uint32_t**) malloc(NUM_DEVICES*sizeof(uint32_t*));

    return NUM_DEVICES;
}

void SendRefWriteRequest (char* seq, size_t start_addr, uint32_t len){

    ref_len = len;
    
    for(int g = 0; g < NUM_DEVICES; g++){

        check_cuda_setDevice(g, "SendRefWriteRequest");

//        const char* d_ref_seq_tmp;
//        check_cuda_malloc((void**)&d_ref_seq_tmp, len*sizeof(char), "tmp_ref_seq"); 

//        check_cuda_memcpy((void*)d_ref_seq_tmp, (void*)(seq + start_addr), len*sizeof(char), cudaMemcpyHostToDevice, "ref_seq");

        check_cuda_malloc((void**)&d_ref_seq[g], len*sizeof(char), "ref_seq"); 

//        compress_string <<<MAX_BLOCKS, MAX_THREADS>>> (len, d_ref_seq_tmp, d_ref_seq[g]);

        char* target_s = (char *) malloc(len*sizeof(char));
        std::memcpy(target_s, seq + start_addr, len);
        pinned_host_vector<int8_t> h_encoded_target(len);
        encode_sequence(h_encoded_target.data(), target_s, len);


        check_cuda_memcpy((void*)d_ref_seq[g], h_encoded_target.data(), len*sizeof(char), cudaMemcpyHostToDevice, "ref_seq");

//        check_cuda_free((void*)d_ref_seq_tmp, "d_ref_seq_tmp");
    }
}

void ClearRef(){

    for(int g = 0; g < NUM_DEVICES; g++){

        check_cuda_setDevice(g, "ClearRef");

        check_cuda_free((void*)d_ref_seq[g], "d_ref_seq");
        check_cuda_free((void*)d_index_table[g], "d_index_table");
        check_cuda_free((void*)d_pos_table[g], "d_pos_table");
    }
}

InitializeInterface_ptr g_InitializeInterface = InitializeInterface;
SendRefWriteRequest_ptr g_SendRefWriteRequest = SendRefWriteRequest;
ClearRef_ptr g_ClearRef = ClearRef;
