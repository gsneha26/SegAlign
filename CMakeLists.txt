cmake_minimum_required (VERSION 3.8)

set(CMAKE_CUDA_COMPILER nvcc)
project(segalign LANGUAGES C CXX CUDA)
find_package(CUDA 9.0 REQUIRED)

set(CMAKE_CUDA_COMPILER_ENV_VAR "CUDACXX")
set(CMAKE_CUDA_HOST_COMPILER_ENV_VAR "CUDAHOSTCXX")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -std=c++11 -l c -lboost_program_options")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_52")

include(${TBB_ROOT}/cmake/TBBBuild.cmake)
tbb_build(TBB_ROOT ${TBB_ROOT} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1)
find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)

find_package(ZLIB REQUIRED)
find_path(ZLIB_INCLUDE_DIR zlib.h)
include_directories(${ZLIB_INCLUDE_DIR})

find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_executable(segalign
    src/DRAM.cpp
    src/GPU.cu
    src/seed_pos_table.cpp
    src/ntcoding.cpp
    src/seeder.cpp
    src/segment_printer.cpp
    src/main.cpp
    )

set_property(TARGET segalign PROPERTY CUDA_STANDARD 11)

# We need to explicitly state that we need all CUDA files in the particle
# library to be built with -dc as the member functions could be called by
# other libraries and executables
set_target_properties(segalign PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_compile_options(segalign PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
                           --generate-line-info>)
target_compile_options(segalign PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)

if(ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
    target_link_libraries (segalign PRIVATE rt stdc++  ${TBB_IMPORTED_TARGETS} ${ZLIB_LIBRARIES} ${CUDA_LIBRARIES})
endif()

if(TARGET zlib::zlib)
    target_link_libraries(segalign PRIVATE rt stdc++  ${TBB_IMPORTED_TARGETS} zlib::zlib ${CUDA_LIBRARIES})
endif()

target_link_libraries(segalign PRIVATE rt stdc++  ${TBB_IMPORTED_TARGETS} ${Boost_LIBRARIES})
