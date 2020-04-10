#include "DRAM.h"
#include <memory>

#include "tbb/scalable_allocator.h"
#include "tbb/tbb.h"

DRAM::DRAM()
	: size(20ull * 1024ull * 1024ull * 1024ull), // 10GB CPU memory
	querySize(0),
	bufferPosition(0)
{

	buffer = (char*)scalable_aligned_malloc(size, 64);
}


DRAM::~DRAM()
{
    scalable_aligned_free(buffer);
}
