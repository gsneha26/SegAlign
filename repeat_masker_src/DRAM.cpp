#include "DRAM.h"
#include <memory>
#include "tbb/scalable_allocator.h"

DRAM::DRAM()
	: size(6ull * 1024ull * 1024ull * 1024ull), // 6GB CPU memory
	seqSize(0),
	bufferPosition(0)
{

	buffer = (char*)scalable_aligned_malloc(size, 64);
}

DRAM::~DRAM(){
    scalable_aligned_free(buffer);
}
