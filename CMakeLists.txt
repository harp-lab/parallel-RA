cmake_minimum_required (VERSION 2.6)

project (BTree)

add_compile_options(--std=c++14 -Wno-strict-aliasing)


find_package(MPI REQUIRED)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -std=c++14")
set (source_dir "${PROJECT_SOURCE_DIR}/src")




file (GLOB source_files_btctp "${source_dir}/hash.cpp" "${source_dir}/vector_buffer.cpp" "${source_dir}/comm.cpp" "${source_dir}/google_btree_relation.cpp" "${source_dir}/all_to_all_comm.cpp" "${source_dir}/parallel_io.cpp" "${source_dir}/parallel_join.cpp" "${source_dir}/parallel_copy.cpp" "${source_dir}/balanced_hash_relation.cpp" "${source_dir}/relation_load_balancer.cpp" "${source_dir}/RA_tasks.cpp" "${source_dir}/lie.cpp" "${source_dir}/transitive_closure_task_parallel.cpp")



add_executable(TC ${source_files_btctp})
INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
TARGET_LINK_LIBRARIES(TC ${MPI_LIBRARIES})


