cmake_minimum_required (VERSION 3.8)

set(CMAKE_CUDA_COMPILER nvcc)
project(wga LANGUAGES C CXX CUDA)
find_package(CUDA 8.0 REQUIRED)
#find_package(Threads REQUIRED)
set(CMAKE_CUDA_COMPILER_ENV_VAR "CUDACXX")
set(CMAKE_CUDA_HOST_COMPILER_ENV_VAR "CUDAHOSTCXX")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -l c -lboost_program_options")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_52")

include(${TBB_ROOT}/cmake/TBBBuild.cmake)
tbb_build(TBB_ROOT ${TBB_ROOT} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1)
find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)


find_package(ZLIB REQUIRED)
FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
# To find and use zlib
find_path(ZLIB_INCLUDE_DIR zlib.h)
include_directories(${ZLIB_INCLUDE_DIR})

FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

add_executable(wga
    DRAM.cpp
    GPU.cu
    seed_pos_table.cpp
    ntcoding.cpp
    seeder.cpp
    segment_printer.cpp
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

TARGET_LINK_LIBRARIES(wga PRIVATE rt stdc++  ${TBB_IMPORTED_TARGETS} pthread ${Boost_LIBRARIES})
