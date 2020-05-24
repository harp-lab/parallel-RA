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

#file (GLOB source_files_btree "${source_dir}/btree.cpp")
#file (GLOB source_files_tc "${source_dir}/transitive_closure.cpp")
#file (GLOB source_files_mpitc "${source_dir}/mpi_transitive_closure.cpp")




#file (GLOB source_files_tom_btc "${source_dir}/RA_tasks.h"  "${source_dir}/parallel_join.h" "${source_dir}/TOM_transitive_closure.cpp" "${source_dir}/TOM_transitive_closure.cpp")
file (GLOB source_files_btctp "${source_dir}/RA_tasks.h" "${source_dir}/lie.h"  "${source_dir}/parallel_join.h" "${source_dir}/balanced_hash_relation.cpp" "${source_dir}/transitive_closure_task_parallel.cpp")
#file (GLOB source_files_btc "${source_dir}/RA_tasks.h" "${source_dir}/lie.h"  "${source_dir}/parallel_join.h" "${source_dir}/balanced_hash_relation.cpp" "${source_dir}/balanced_transitive_closure.cpp")
#file (GLOB source_files_btcdp "${source_dir}/RA_tasks.h"  "${source_dir}/parallel_join.h" "${source_dir}/balanced_hash_relation.cpp" "${source_dir}/balanced_transitive_closure_data_parallel.cpp")
#file (GLOB source_files_btcdp2 "${source_dir}/RA_tasks.h"  "${source_dir}/parallel_join.h" "${source_dir}/balanced_hash_relation.cpp" "${source_dir}/balanced_transitive_closure_data_parallel2.cpp")

#add_library(hash_relation "${source_dir}/balanced_hash_relation.cpp")



#file (GLOB source_files_btc_tjt "${source_dir}/balanced_transitive_closure_t_join_t.cpp")
#file (GLOB source_files_btc_data_parallel "${source_dir}/balanced_transitive_closure_data_parallel.cpp")
#file (GLOB source_files_copy "${source_dir}/balanced_copy.cpp")
#file (GLOB source_files_union "${source_dir}/union.cpp")
#file (GLOB source_files_aggregated_union "${source_dir}/aggregated_union.cpp")
#file (GLOB source_files_new_join "${source_dir}/new_join.cpp")
#file (GLOB source_files_join "${source_dir}/join.cpp")
#file (GLOB source_files_ata "${source_dir}/all_to_all.cpp")
#file (GLOB source_files_mpitc_map "${source_dir}/mpi_tc_unordered_map.cpp")
#file (GLOB source_files_mpitc_google_map "${source_dir}/mpi_tc_google_map.cpp")
#file (GLOB source_files_mpi_double_hash_tc_google_map "${source_dir}/mpi_double_hash_tc_google_map.cpp")
#file (GLOB source_files_mpi_buffered_double_hash_tc_google_map "${source_dir}/mpi_buffered_double_hash_tc_google_map.cpp")
#file (GLOB source_files_balanced_hash_tree_TC "${source_dir}/balanced_hash_tree_TC.cpp")
#file (GLOB source_files_mpi_load_balancer_test "${source_dir}/mpi_load_balancer_test.cpp")
#file (GLOB source_files_double_hash_tc_google_map "${source_dir}/double_hash_tc_google_map.cpp")
#file (GLOB source_files_otc "${source_dir}/optimized_transitive_closure.cpp")
#file (GLOB source_files_relation "${source_dir}/relation.cpp")
#file (GLOB source_files_dup_relation "${source_dir}/dup_relation.cpp")

#add_executable(BTree ${source_files_btree})

#add_executable(TC ${source_files_tc})

#add_executable(ata ${source_files_ata})
#INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
#TARGET_LINK_LIBRARIES(ata ${MPI_LIBRARIES})

#add_executable(mpiTC ${source_files_mpitc})
#INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
#TARGET_LINK_LIBRARIES(mpiTC ${MPI_LIBRARIES})



#add_executable(tom_bTC ${source_files_tom_btc})
#INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
#TARGET_LINK_LIBRARIES(tom_bTC ${MPI_LIBRARIES})


add_executable(TC ${source_files_btctp})
INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
TARGET_LINK_LIBRARIES(TC ${MPI_LIBRARIES})

#add_executable(bTC ${source_files_btc})
#INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
#TARGET_LINK_LIBRARIES(bTC ${MPI_LIBRARIES})

#add_executable(bTCdp ${source_files_btcdp})
#INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
#TARGET_LINK_LIBRARIES(bTCdp ${MPI_LIBRARIES})

#add_executable(bTCdp2 ${source_files_btcdp2})
#INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
#TARGET_LINK_LIBRARIES(bTCdp2 ${MPI_LIBRARIES})





#add_executable(b_copy ${source_files_copy})
#INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
#TARGET_LINK_LIBRARIES(b_copy ${MPI_LIBRARIES})

#add_executable(bTC_tjt ${source_files_btc_tjt})
#INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
#TARGET_LINK_LIBRARIES(bTC_tjt ${MPI_LIBRARIES})

#add_executable(bTC_dp ${source_files_btc_data_parallel})
#INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
#TARGET_LINK_LIBRARIES(bTC_dp ${MPI_LIBRARIES})


#add_executable(union ${source_files_union})
#INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
#TARGET_LINK_LIBRARIES(union ${MPI_LIBRARIES})


#add_executable(aggregated_union ${source_files_aggregated_union})
#INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
#TARGET_LINK_LIBRARIES(aggregated_union ${MPI_LIBRARIES})


#add_executable(join ${source_files_join})
#INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
#TARGET_LINK_LIBRARIES(join ${MPI_LIBRARIES})

#add_executable(new_join ${source_files_new_join})
#INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
#TARGET_LINK_LIBRARIES(new_join ${MPI_LIBRARIES})


#add_executable(mpiTCmap ${source_files_mpitc_map})
#INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
#TARGET_LINK_LIBRARIES(mpiTCmap ${MPI_LIBRARIES})


#add_executable(mpiTCGoogleMap ${source_files_mpitc_google_map})
#INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
#TARGET_LINK_LIBRARIES(mpiTCGoogleMap ${MPI_LIBRARIES})

#add_executable(mpiTC_double_hash_googleMap ${source_files_mpi_double_hash_tc_google_map})
#INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
#TARGET_LINK_LIBRARIES(mpiTC_double_hash_googleMap ${MPI_LIBRARIES})


#add_executable(mpiTC_buffered_double_hash_googleMap ${source_files_mpi_buffered_double_hash_tc_google_map})
#INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
#TARGET_LINK_LIBRARIES(mpiTC_buffered_double_hash_googleMap ${MPI_LIBRARIES})

#add_executable(balanced_hash_tree_TC ${source_files_balanced_hash_tree_TC})
#INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
#TARGET_LINK_LIBRARIES(balanced_hash_tree_TC ${MPI_LIBRARIES})


#add_executable(mpi_load_balancer_test ${source_files_mpi_load_balancer_test})
#INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
#TARGET_LINK_LIBRARIES(mpi_load_balancer_test ${MPI_LIBRARIES})

#add_executable(TC_double_hash_googleMap ${source_files_double_hash_tc_google_map})

#add_executable(oTC ${source_files_otc})

#add_executable(relation ${source_files_relation})
#add_executable(dup_relation ${source_files_dup_relation})
