cmake_minimum_required (VERSION 3.8)

project(segalign LANGUAGES CXX CUDA)
project(segalign_repeat_masker LANGUAGES C CXX CUDA)

set(CMAKE_CUDA_COMPILER nvcc)
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
    src/seed_filter.cu
    src/ntcoding.cpp
    src/seeder.cpp
    src/segment_printer.cpp
    src/main.cpp
    )

add_executable(segalign_repeat_masker
    repeat_masker_src/DRAM.cpp
    repeat_masker_src/seed_filter.cu
    repeat_masker_src/ntcoding.cpp
    repeat_masker_src/seeder.cpp
    repeat_masker_src/segment_printer.cpp
    repeat_masker_src/main.cpp
    )

set(LIBS rt stdc++ ${TBB_IMPORTED_TARGETS} ${Boost_LIBRARIES} ${CUDA_LIBRARIES} ${ZLIB_LIBRARIES})

set_property(TARGET segalign PROPERTY CUDA_STANDARD 11)

set_target_properties(segalign PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_compile_options(segalign PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)

target_link_libraries(segalign PRIVATE ${LIBS})

set_property(TARGET segalign_repeat_masker PROPERTY CUDA_STANDARD 11)

set_target_properties(segalign_repeat_masker PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_compile_options(segalign_repeat_masker PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)

target_link_libraries(segalign_repeat_masker PRIVATE ${LIBS})

#if(ZLIB_FOUND)
#    target_link_libraries (segalign PRIVATE ${ZLIB_LIBRARIES})
#    target_link_libraries (segalign_repeat_masker PRIVATE ${ZLIB_LIBRARIES})
#endif()
#
#if(TARGET zlib::zlib)
#    target_link_libraries(segalign PRIVATE zlib::zlib )
#    target_link_libraries(segalign_repeat_masker PRIVATE zlib::zlib )
#endif()
#
#message(${ZLIB_LIBRARIES})
#message(${Boost_LIBRARIES})
#message(${CUDA_LIBRARIES})
#message(${TBB_IMPORTED_TARGETS})
#message(zlib::zlib)
#message(${LIBS})
