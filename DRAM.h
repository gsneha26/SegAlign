#pragma once
#include <cstddef>

#define WORD_SIZE 128

class DRAM
{
public:
	char* buffer;
	std::size_t size;

	std::size_t referenceSize;
	std::size_t bufferPosition;
public:
	DRAM();
	~DRAM();
};

