cmake_minimum_required (VERSION 3.8)
cmake_policy(SET CMP0104 NEW)
cmake_policy(SET CMP0128 NEW)
cmake_policy(SET CMP0146 NEW)

execute_process(COMMAND bash "-c" "${CMAKE_SOURCE_DIR}/scripts/get-cuda-arches.bash" OUTPUT_VARIABLE SEGALIGN_CUDA_ARCHITECTURES)
set(CMAKE_CUDA_ARCHITECTURES "${SEGALIGN_CUDA_ARCHITECTURES}")

project(segalign LANGUAGES C CXX CUDA)
project(segalign_repeat_masker LANGUAGES C CXX CUDA)

set(CMAKE_CUDA_COMPILER nvcc)

set(CMAKE_CUDA_COMPILER_ENV_VAR "CUDACXX")
set(CMAKE_CUDA_HOST_COMPILER_ENV_VAR "CUDAHOSTCXX")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -std=c++14 -l c -lboost_program_options")

find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)

find_package(ZLIB REQUIRED)
find_path(ZLIB_INCLUDE_DIR zlib.h)
include_directories(${ZLIB_INCLUDE_DIR})

find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_executable(segalign
    common/DRAM.cpp
    common/ntcoding.cpp
    common/seed_filter_interface.cu
    common/seed_pos_table.cu
    common/utilities.c
    common/dna_utilities.c
    common/scoring.c
    src/seed_filter.cu
    src/seeder.cpp
    src/segment_printer.cpp
    src/main.cpp
    )

add_executable(segalign_repeat_masker
    common/DRAM.cpp
    common/ntcoding.cpp
    common/seed_filter_interface.cu
    common/seed_pos_table.cu
    common/utilities.c
    common/dna_utilities.c
    common/scoring.c
    repeat_masker_src/seed_filter.cu
    repeat_masker_src/seeder.cpp
    repeat_masker_src/segment_printer.cpp
    repeat_masker_src/main.cpp
    )

set(LIBS rt stdc++ ${TBB_IMPORTED_TARGETS} ${Boost_LIBRARIES} ${CUDA_LIBRARIES} ${ZLIB_LIBRARIES})

set_target_properties(
    segalign 
    segalign_repeat_masker
    PROPERTIES 
      CUDA_SEPARABLE_COMPILATION ON
      CUDA_STANDARD 14
      )

target_compile_options(segalign PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
target_compile_options(segalign_repeat_masker PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)

target_link_libraries(segalign PRIVATE ${LIBS})
target_link_libraries(segalign_repeat_masker PRIVATE ${LIBS})

target_include_directories(segalign PRIVATE ${PROJECT_SOURCE_DIR}/common/ ${PROJECT_SOURCE_DIR}/src/)
target_include_directories(segalign_repeat_masker PRIVATE ${PROJECT_SOURCE_DIR}/common/ ${PROJECT_SOURCE_DIR}/repeat_masker_src/)
