#include <stdio.h>

// wrap of cudaSetDevice error checking in one place.  
static inline void check_cuda_setDevice(int device_id, const char* tag) {
    cudaError_t err = cudaSetDevice(device_id);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaSetDevice failed for device %d in %s failed with error \" %s \" \n", device_id, tag, cudaGetErrorString(err));
        exit(11);
    }
}

// wrap of cudaMalloc error checking in one place.  
static inline void check_cuda_malloc(void** buf, size_t bytes, const char* tag) {
    cudaError_t err = cudaMalloc(buf, bytes);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMalloc of %lu bytes for %s failed with error \" %s \" \n", bytes, tag, cudaGetErrorString(err));
        exit(12);
    }
}
	 
// wrap of cudaMemcpy error checking in one place.  
static inline void check_cuda_memcpy(void* dst_buf, void* src_buf, size_t bytes, cudaMemcpyKind kind, const char* tag) {
    cudaError_t err = cudaMemcpy(dst_buf, src_buf, bytes, kind);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy of %lu bytes for %s failed with error \" %s \" \n", bytes, tag, cudaGetErrorString(err));
        exit(13);
    }
}
	 
// wrap of cudaFree error checking in one place.  
static inline void check_cuda_free(void* buf, const char* tag) {
    cudaError_t err = cudaFree(buf);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaFree for %s failed with error \" %s \" \n", tag, cudaGetErrorString(err));
        exit(14);
    }
}
