#ifndef RAM_TASKS_H
#define RAM_TASKS_H

#include <mpi.h>
#include <vector>
#include "balanced_hash_relation.h"
#include "balanced_parallel_io.h"
#include "parallel_join.h"
//#include "google_btree_relation.h"
#include "btree/btree_map.h"
#include <queue>

class RAM
{
private:
    std::vector<parallel_RA> RA_list;

    u64 *clique_buf_output_size;
    u64 **clique_buf_output;

    std::queue<u64> *clique_output_queue;

    vector_buffer** local_join_output;
    int **local_join_output_size;

    u64** outer_hash_data;
    int *outer_hash_buffer_size;

    mpi_comm mcomm;


public:

    RAM(mpi_comm mcomm)
    {
        this->mcomm = mcomm;
    }

    void push_back(parallel_RA pj)
    {
        RA_list.push_back(pj);
    }


    // Print relation size
    // this function is hacky right now
    void print_full()
    {
        for (std::vector<parallel_RA>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
        {
            parallel_RA current_ra = *it;
            if (current_ra.get_RA_type() == COPY)
            {
                relation* input = current_ra.get_copy_input();
                relation* output = current_ra.get_copy_output();

                std::cout << "Input 0 " << current_ra.full_count(input) << std::endl;
                std::cout << "Input 1 " << current_ra.full_count(output) << std::endl;
            }
        }
        return;
    }


    // check for fixed point
    bool check_for_fixed_point()
    {
        bool fixed_point = true;
        for (std::vector<parallel_RA>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
        {
            parallel_RA current_ra = *it;
            if (current_ra.get_RA_type() == COPY)
            {
                relation* output = current_ra.get_copy_output();
                fixed_point = fixed_point &  current_ra.fixed_point_check(output);
            }
            else if (current_ra.get_RA_type() == JOIN)
            {
                relation* output = current_ra.get_join_output();
                fixed_point = fixed_point & current_ra.fixed_point_check(output);
            }
        }

        return fixed_point;
    }


    // Clique Comm
    u64 clique_comm()
    {
        u64 total_data_moved = 0;
        u32 counter = 0;
        u32 RA_count = RA_list.size();

        clique_buf_output_size = new u64[RA_count];
        clique_buf_output = new u64*[RA_count];
        clique_output_queue = new std::queue<u64>[RA_count];

        for (std::vector<parallel_RA>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
        {
            parallel_RA current_ra = *it;
            if (current_ra.get_RA_type() == COPY)
            {
                relation* input = current_ra.get_copy_input();
                relation* output = current_ra.get_copy_output();

                if (current_ra.get_copy_input0_graph_type() == FULL)
                {
                    current_ra.full_clique_comm(input, output, &clique_buf_output_size[counter], &clique_buf_output[counter]);
                    total_data_moved = total_data_moved + clique_buf_output_size[counter];
                }
                else if (current_ra.get_copy_input0_graph_type() == DELTA)
                {
                    current_ra.delta_clique_comm(input, output, &clique_buf_output_size[counter], &clique_buf_output[counter]);
                    total_data_moved = total_data_moved + clique_buf_output_size[counter];
                }
            }

            else if (current_ra.get_RA_type() == JOIN)
            {
                relation* input = current_ra.get_join_input1();
                relation* output = current_ra.get_join_input0();

                if (current_ra.get_join_input1_graph_type() == FULL)
                {
                    current_ra.full_clique_comm(input, output, &clique_buf_output_size[counter], &clique_buf_output[counter]);
                    total_data_moved = total_data_moved + clique_buf_output_size[counter];
                }
                else if (current_ra.get_join_input1_graph_type() == DELTA)
                {
                    current_ra.delta_clique_comm(input, output, &clique_buf_output_size[counter], &clique_buf_output[counter]);
                    total_data_moved = total_data_moved + clique_buf_output_size[counter];
                }
            }

            for (u64 i = 0; i < clique_buf_output_size[counter]; i++)
                clique_output_queue[counter].push(clique_buf_output[counter][i]);

            counter++;
        }



        return total_data_moved;
    }

    u64 local_join()
    {
        u32 total_join_tuples = 0;
        u32 RA_count = RA_list.size();
        u32 nprocs = (u32)mcomm.get_nprocs();

        local_join_output = new vector_buffer*[RA_count];
        local_join_output_size = new int*[RA_count];
        for (u32 i = 0; i < RA_count; i++)
        {
            local_join_output[i] = new vector_buffer[nprocs];
            for (u32 j = 0; j < nprocs; ++j) {
                local_join_output[i][j] = vector_buffer_create_empty();
            }

            local_join_output_size[i] = new int[nprocs];
            memset(local_join_output_size[i], 0, nprocs * sizeof(int));
        }

        u32 counter = 0;
        for (std::vector<parallel_RA>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
        {
            parallel_RA current_ra = *it;
            if (current_ra.get_RA_type() == COPY)
            {
                relation* output = current_ra.get_copy_output();
                current_ra.local_copy(clique_buf_output_size[counter], 2, clique_buf_output[counter], output, local_join_output[counter], local_join_output_size[counter]);
            }
            else if (current_ra.get_RA_type() == JOIN)
            {
                relation* input0 = current_ra.get_join_input0();
                relation* output = current_ra.get_join_output();

                if (current_ra.get_join_input0_graph_type() == FULL)
                {
                    total_join_tuples = total_join_tuples + current_ra.local_join_full(clique_buf_output_size[counter], 2, clique_buf_output[counter], input0, output, local_join_output[counter], local_join_output_size[counter], current_ra.get_projection_index());
                }
                else if (current_ra.get_join_input0_graph_type() == DELTA)
                {
                    total_join_tuples = total_join_tuples + current_ra.local_join_delta(clique_buf_output_size[counter], 2, clique_buf_output[counter], input0, output, local_join_output[counter], local_join_output_size[counter], current_ra.get_projection_index());
                }
            }
            counter++;
        }

        delete[] clique_buf_output_size;
        delete[] clique_buf_output;


        return total_join_tuples;
    }


    int all_to_all()
    {
        u64 total_all_to_all = 0;

        u32 RA_count = RA_list.size();

        outer_hash_data = new u64*[RA_count];
        outer_hash_buffer_size = new int[RA_count];

        u32 counter = 0;
        for (std::vector<parallel_RA>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
        {
            parallel_RA current_ra = *it;
            current_ra.all_to_all(local_join_output[counter], local_join_output_size[counter], &(outer_hash_buffer_size[counter]), &(outer_hash_data[counter]));
            total_all_to_all = total_all_to_all + outer_hash_buffer_size[counter];
            counter++;
        }

        for (u32 i = 0; i < RA_count; i++)
        {
            delete[] local_join_output[i];
            delete[] local_join_output_size[i];
        }
        delete[] local_join_output;
        delete[] local_join_output_size;


        return total_all_to_all;
    }



    u64 local_insert_in_full()
    {
        u64 local_insert_in_full = 0;
        for (std::vector<parallel_RA>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
        {
            parallel_RA current_ra = *it;
            if (current_ra.get_RA_type() == COPY)
            {
                relation* output = current_ra.get_copy_output();
                local_insert_in_full = local_insert_in_full + current_ra.local_inserts_in_full(output);
            }
            else if (current_ra.get_RA_type() == JOIN)
            {
                relation* output = current_ra.get_join_output();
                local_insert_in_full = local_insert_in_full + current_ra.local_inserts_in_full(output);
            }
        }

        return local_insert_in_full;
    }



    u64 local_insert_in_delta()
    {
        u64 local_insert_in_delta = 0;
        u32 counter = 0;
        for (std::vector<parallel_RA>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
        {
            parallel_RA current_ra = *it;
            if (current_ra.get_RA_type() == COPY)
            {
                relation* output = current_ra.get_copy_output();
                local_insert_in_delta = local_insert_in_delta + current_ra.local_insert_in_delta(output, outer_hash_buffer_size[counter], outer_hash_data[counter]);
            }
            else if (current_ra.get_RA_type() == JOIN)
            {
                relation* output = current_ra.get_join_output();
                local_insert_in_delta = local_insert_in_delta + current_ra.local_insert_in_delta(output, outer_hash_buffer_size[counter], outer_hash_data[counter]);
            }
            counter++;
        }

        delete[] outer_hash_data;
        delete[] outer_hash_buffer_size;

        return local_insert_in_delta;
    }


    void execute()
    {
        u32 iteration = 0;
        u32 clique_comm_count = 0;
        u32 local_join_count = 0;
        u32 all_to_all_count = 0;
        u32 insert_in_full_count = 0;
        u32 insert_in_delta = 0;

        u32 running_clique_comm_count = 0;
        u32 running_local_join_count = 0;
        u32 running_all_to_all_count = 0;
        u32 running_insert_in_full_count = 0;
        u32 running_insert_in_delta = 0;

        double clique_start = 0, clique_end = 0;
        double local_join_start = 0, local_join_end = 0;
        double all_to_all_start = 0, all_to_all_end = 0;
        double insert_full_start = 0, insert_full_end = 0;
        double insert_delta_start = 0, insert_delta_end = 0;
        double verify_start = 0, verify_end = 0;

        double iteration_time = 0;
        double running_time = 0;

        while (true)
        {
            clique_start = MPI_Wtime();
            clique_comm_count = clique_comm();
            running_clique_comm_count = running_clique_comm_count + clique_comm_count;
            clique_end = MPI_Wtime();

            local_join_start = MPI_Wtime();
            local_join_count = local_join();
            running_local_join_count = running_local_join_count + local_join_count;
            local_join_end = MPI_Wtime();

            std::cout << "All to All" << std::endl;

            all_to_all_start = MPI_Wtime();
            all_to_all_count = all_to_all();
            running_all_to_all_count = running_all_to_all_count + all_to_all_count;
            all_to_all_end = MPI_Wtime();

            insert_full_start = MPI_Wtime();
            insert_in_full_count = local_insert_in_full();
            running_insert_in_full_count = running_insert_in_full_count + insert_in_full_count;
            insert_full_end = MPI_Wtime();

            insert_delta_start = MPI_Wtime();
            insert_in_delta = local_insert_in_delta();
            running_insert_in_delta = running_insert_in_delta + insert_in_delta;
            insert_delta_end = MPI_Wtime();


            //std::cout << "ITERATION [" << iteration <<"] " << clique_comm_count << " " << local_join_count << " " << all_to_all_count << " " << insert_in_full_count << " " << insert_in_delta << std::endl;
            //std::cout << "RUNNING ITERATION [" << iteration <<"] " << running_clique_comm_count << " " << running_local_join_count << " " << running_all_to_all_count << " " << running_insert_in_full_count << " " << running_insert_in_delta << std::endl;


            verify_start = MPI_Wtime();
            if (check_for_fixed_point() == true)
            {
                verify_end = MPI_Wtime();
                iteration_time = (verify_end - verify_start) + (insert_delta_end - insert_delta_start) + (insert_full_end - insert_full_start) + (all_to_all_end - all_to_all_start) + (local_join_end - local_join_start) + (clique_end - clique_start);
                running_time = running_time + iteration_time;

                std::cout << " T [" << iteration << "] " << running_time << " " << iteration_time
                          << " Clique " <<  (clique_end - clique_start)
                          << " Local Join " <<  (local_join_end - local_join_start)
                          << " All to All " <<  (all_to_all_end - all_to_all_start)
                          << " Insert Full " <<  (insert_full_end - insert_full_start)
                          << " Insert delta " <<  (insert_delta_end - insert_delta_start)
                          << " Verify " <<  (verify_end - verify_start)
                          << std::endl
                          << std::endl;


                break;
            }
            verify_end = MPI_Wtime();

            iteration_time = (verify_end - verify_start) + (insert_delta_end - insert_delta_start) + (insert_full_end - insert_full_start) + (all_to_all_end - all_to_all_start) + (local_join_end - local_join_start) + (clique_end - clique_start);
            running_time = running_time + iteration_time;

            std::cout << " F [" << iteration << "] " << running_time << " " << iteration_time
                      << " Clique " <<  (clique_end - clique_start)
                      << " Local Join " <<  (local_join_end - local_join_start)
                      << " All to All " <<  (all_to_all_end - all_to_all_start)
                      << " Insert Full " <<  (insert_full_end - insert_full_start)
                      << " Insert delta " <<  (insert_delta_end - insert_delta_start)
                      << " Verify " <<  (verify_end - verify_start)
                      << std::endl
                      << std::endl;

            iteration++;
        }

        print_full();
    }

};

#endif
