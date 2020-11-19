cmake_minimum_required (VERSION 3.9)

project (BTree)

add_compile_options(--std=c++17 -Wno-strict-aliasing)

find_package(MPI REQUIRED)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -std=c++14")
set (source_dir "${PROJECT_SOURCE_DIR}/src")
set (tests_dir "${PROJECT_SOURCE_DIR}/tests")
set (utility_dir "${PROJECT_SOURCE_DIR}/utility")

file (GLOB source_files_parallel_RA "${source_dir}/parallel_RA_inc.h" "${source_dir}/log/logger.cpp" "${source_dir}/hash/hash.h" "${source_dir}/buffer/vector_buffer.cpp" "${source_dir}/comm/comm.cpp" "${source_dir}/relation/google_btree_relation.cpp" "${source_dir}/comm/all_to_all_comm.cpp" "${source_dir}/comm/all_to_allv_comm.cpp" "${source_dir}/IO/parallel_io.cpp" "${source_dir}/RA/parallel_join.cpp" "${source_dir}/comm/intra_bucket_comm.cpp" "${source_dir}/RA/parallel_copy.cpp" "${source_dir}/RA/parallel_copy_filter.cpp" "${source_dir}/RA/parallel_RA.h" "${source_dir}/RA/parallel_acopy.cpp" "${source_dir}/relation/balanced_hash_relation.cpp" "${source_dir}/relation/relation_load_balancer.cpp" "${source_dir}/RAM/RA_tasks.cpp" "${source_dir}/lie/lie.cpp")
file (GLOB source_files_tc "${tests_dir}/transitive_closure.cpp")
file (GLOB source_files_fv "${tests_dir}/freevars.cpp")
file (GLOB source_files_kcfa "${tests_dir}/kcfa-tiny.cpp")

ADD_LIBRARY(parallel_RA "${source_files_parallel_RA}")

add_executable(TC ${source_files_tc})
INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
TARGET_LINK_LIBRARIES(TC ${MPI_LIBRARIES} parallel_RA)

add_executable(freevars ${source_files_fv})
INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
TARGET_LINK_LIBRARIES(freevars ${MPI_LIBRARIES} parallel_RA)

add_executable(kcfa ${source_files_kcfa})
INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
TARGET_LINK_LIBRARIES(kcfa ${MPI_LIBRARIES} parallel_RA)

add_executable(bin_parser "${utility_dir}/binary_parser.cpp")
add_executable(tsv_to_bin "${utility_dir}/tsv_to_bin.cpp")
