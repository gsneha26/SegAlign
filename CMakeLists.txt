cmake_minimum_required (VERSION 3.8)                                                                                                                                                     

#project(wga)

set(CMAKE_CUDA_COMPILER /usr/local/cuda-10.2/bin/nvcc)
project(wga LANGUAGES C CXX CUDA)
find_package(CUDA 8.0 REQUIRED)
#find_package(Threads REQUIRED)
set(CMAKE_CUDA_COMPILER_ENV_VAR "CUDACXX")
set(CMAKE_CUDA_HOST_COMPILER_ENV_VAR "CUDAHOSTCXX")

#set(CMAKE_CXX_COMPILER "nvcc")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -g -pg -std=c++11")
#set(CMAKE_NVCC_FLAGS "${CMAKE_NVCC_FLAGS} -std=c++11 –arch=sm_60")
#set(CMAKE_NVCC_FLAGS "${CMAKE_NVCC_FLAGS} -std=c++11 -g -G -lineinfo –arch=sm_60")

include(${TBB_ROOT}/cmake/TBBBuild.cmake)
tbb_build(TBB_ROOT ${TBB_ROOT} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1)
find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)

#set (XILINX_LINK_LIBS libxilinxopencl.so)
#link_directories(${XILINX_SDX}/runtime/lib/x86_64 ${LD_LIBRARY_PATH})

find_package(ZLIB REQUIRED)
# To find and use zlib
find_path(ZLIB_INCLUDE_DIR zlib.h)
include_directories(${ZLIB_INCLUDE_DIR})

add_executable(wga
    Chameleon.cpp
    ConfigFile.cpp
    DRAM.cpp
    GPU.cu
    seed_pos_table.cpp
    ntcoding.cpp
    seeder.cpp
    filter.cpp
    extender.cpp
    maf_printer.cpp
    main.cpp
    )
set_property(TARGET wga PROPERTY CUDA_STANDARD 11)

#add_executable(wga main.cpp)

# We need to explicitly state that we need all CUDA files in the particle
# library to be built with -dc as the member functions could be called by
# other libraries and executables
set_target_properties(wga PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_compile_options(wga PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
                           --generate-line-info>)


if(ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
    target_link_libraries (wga PRIVATE rt stdc++  ${TBB_IMPORTED_TARGETS} pthread ${ZLIB_LIBRARIES} ${CUDA_LIBRARIES})
endif()

if(TARGET zlib::zlib)
    target_link_libraries(wga PRIVATE rt stdc++  ${TBB_IMPORTED_TARGETS} pthread zlib::zlib ${CUDA_LIBRARIES})
endif()

#execute_process (
#    COMMAND bash -c "${XILINX_SDX}/bin/emconfigutil --od . --nd 1  --platform ${AWS_PLATFORM}")     

