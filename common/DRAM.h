#pragma once
#include <cstddef>

class DRAM{
    public:
        char* buffer;
        std::size_t size;
        std::size_t seqSize;
        std::size_t bufferPosition;
        DRAM();
        ~DRAM();
};
