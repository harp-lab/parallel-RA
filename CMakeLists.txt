cmake_minimum_required (VERSION 2.6)
project (BTree)

find_package(MPI REQUIRED)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -std=c++14")
set (source_dir "${PROJECT_SOURCE_DIR}/src")

file (GLOB source_files_btree "${source_dir}/btree.cpp")
file (GLOB source_files_tc "${source_dir}/transitive_closure.cpp")
file (GLOB source_files_mpitc "${source_dir}/mpi_transitive_closure.cpp")
file (GLOB source_files_otc "${source_dir}/optimized_transitive_closure.cpp")
file (GLOB source_files_relation "${source_dir}/relation.cpp")

add_executable(BTree ${source_files_btree})

add_executable(TC ${source_files_tc})

add_executable(mpiTC ${source_files_mpitc})
INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
TARGET_LINK_LIBRARIES(mpiTC ${MPI_LIBRARIES})

#add_executable(oTC ${source_files_otc})

add_executable(relation ${source_files_relation})
