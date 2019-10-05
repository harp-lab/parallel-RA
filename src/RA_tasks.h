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

u32 running_total = 0;
bool threshold_enabled = false;
u32 threshold = 10;

class RAM
{
private:
    std::vector<parallel_RA*> RA_list;

    u64 *clique_buf_output_size;
    u64 **clique_buf_output;

    vector_buffer** local_join_output;
    int **local_join_output_size;

    mpi_comm mcomm;


public:

    RAM(mpi_comm mcomm)
    {
        this->mcomm = mcomm;
    }

    void push_back(parallel_RA* pj)
    {
        RA_list.push_back(pj);
    }


    // Print relation size
    // this function is hacky right now
    void print_full()
    {
        int rank = mcomm.get_rank();

        for (std::vector<parallel_RA*>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
        {
            parallel_RA* current_ra = *it;
            if (current_ra->get_RA_type() == COPY)
            {
                relation* input = current_ra->get_copy_input();
                relation* output = current_ra->get_copy_output();
                u32 s0 = current_ra->full_count(input);
                u32 s1 = current_ra->full_count(output);

                if (rank == 0)
                {
                std::cout << "Input 0 " << s0 << std::endl;
                std::cout << "Input 1 " << s1 << std::endl;
                }
            }
            else if (current_ra->get_RA_type() == JOIN)
            {
                relation* output = current_ra->get_join_output();
                int s0 = current_ra->full_count(output);

                if (rank == 0)
                {
                    std::cout << "Output 1 " << s0 << std::endl;
                }
            }
        }
        return;
    }


    // check for fixed point
    bool check_for_fixed_point(bool local_join_status)
    {
        if (local_join_status == false)
            return false;

        bool fixed_point = true;
        u32 counter = 0;
        for (std::vector<parallel_RA*>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
        {
            parallel_RA* current_ra = *it;
            if (current_ra->get_RA_type() == COPY)
            {
                relation* output = current_ra->get_copy_output();
                fixed_point = fixed_point &  current_ra->fixed_point_check(output);
            }
            else if (current_ra->get_RA_type() == JOIN)
            {
                relation* input0 = current_ra->get_join_input0();
                relation* input1 = current_ra->get_join_input1();
                relation* output = current_ra->get_join_output();
                fixed_point = fixed_point & current_ra->fixed_point_check(input0)
                                          & current_ra->fixed_point_check(input1)
                                          & current_ra->fixed_point_check(output);

                //if (threshold_enabled == true)
                //{
                    //std::cout << "VVVVVVVVVVV " << clique_output_queue[counter]->size() << std::endl;
                //    if (clique_output_queue[counter]->size() != 0)
                //        return false;
                //}
            }
            counter++;
        }

        return fixed_point;
    }


    // Clique Comm
    u64 clique_comm(bool local_join_status)
    {
        if (local_join_status == false)
            return 0;

        u64 total_data_moved = 0;
        u32 counter = 0;
        u32 RA_count = RA_list.size();

        clique_buf_output_size = new u64[RA_count];
        clique_buf_output = new u64*[RA_count];


        for (std::vector<parallel_RA*>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
        {
            parallel_RA* current_ra = *it;
            if (current_ra->get_RA_type() == COPY)
                continue;
            else if (current_ra->get_RA_type() == JOIN)
            {

                if (current_ra->get_join_input0_graph_type() == DELTA && current_ra->get_join_input1_graph_type() == DELTA)
                {
                    // room for optimization
                    relation* input = current_ra->get_join_input1();
                    relation* output = current_ra->get_join_input0();

                    current_ra->clique_comm(input->get_arity(), input->get_delta(), input->get_delta_element_count(), input->get_distinct_sub_bucket_rank_count(), input->get_distinct_sub_bucket_rank(), input->get_bucket_map(), output, &clique_buf_output_size[counter], &clique_buf_output[counter]);
                    total_data_moved = total_data_moved + clique_buf_output_size[counter];
                }
                else if (current_ra->get_join_input0_graph_type() == DELTA && current_ra->get_join_input1_graph_type() == FULL)
                {
                    relation* input = current_ra->get_join_input0();
                    relation* output = current_ra->get_join_input1();

                    current_ra->clique_comm(input->get_arity(), input->get_delta(), input->get_delta_element_count(), input->get_distinct_sub_bucket_rank_count(), input->get_distinct_sub_bucket_rank(), input->get_bucket_map(), output, &clique_buf_output_size[counter], &clique_buf_output[counter]);
                    total_data_moved = total_data_moved + clique_buf_output_size[counter];

                }
                else if (current_ra->get_join_input0_graph_type() == FULL && current_ra->get_join_input1_graph_type() == DELTA)
                {
                    relation* input = current_ra->get_join_input1();
                    relation* output = current_ra->get_join_input0();

                    current_ra->clique_comm(input->get_arity(), input->get_delta(), input->get_delta_element_count(), input->get_distinct_sub_bucket_rank_count(), input->get_distinct_sub_bucket_rank(), input->get_bucket_map(), output, &clique_buf_output_size[counter], &clique_buf_output[counter]);
                    total_data_moved = total_data_moved + clique_buf_output_size[counter];
                }

                //for (u64 i = 0; i < clique_buf_output_size[counter]; i++)
                //    clique_output_queue[counter]->push(clique_buf_output[counter][i]);
                //std::cout << "[Queue] pushing " << clique_buf_output_size[counter]/2 << std::endl;
            }
            counter++;
        }

        return total_data_moved;
    }



    bool local_join(u32* local_join_count, int* offset, int iteration)
    {
        u32 join_tuples = 0;
        u32 total_join_tuples = 0;
        u32 RA_count = RA_list.size();
        u32 nprocs = (u32)mcomm.get_nprocs();

        bool join_completed = true;

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
        for (std::vector<parallel_RA*>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
        {
            parallel_RA* current_ra = *it;
            if (current_ra->get_RA_type() == COPY)
            {
                relation* output = current_ra->get_copy_output();
                relation* input = current_ra->get_copy_input();
                current_ra->local_copy(input, output, local_join_output[counter], local_join_output_size[counter]);
            }
            else if (current_ra->get_RA_type() == JOIN)
            {
                    //std::cout << "OFFSET " << offset[counter] << std::endl;
                    if (current_ra->get_join_input0_graph_type() == DELTA && current_ra->get_join_input1_graph_type() == DELTA)
                    {
                        relation* input0 = current_ra->get_join_input0();
                        relation* output = current_ra->get_join_output();

                        join_completed = join_completed & current_ra->local_join(threshold, clique_buf_output_size[counter], &(offset[counter]), 2, clique_buf_output[counter], input0->get_delta(), input0->get_delta_element_count(), output, local_join_output[counter], local_join_output_size[counter], 0, 1, &join_tuples, iteration);
                        total_join_tuples = total_join_tuples + join_tuples;
                    }
                    else if (current_ra->get_join_input0_graph_type() == DELTA && current_ra->get_join_input1_graph_type() == FULL)
                    {
                        relation* input1 = current_ra->get_join_input1();
                        relation* output = current_ra->get_join_output();

                        join_completed = join_completed & current_ra->local_join(threshold, clique_buf_output_size[counter], &(offset[counter]), 2, clique_buf_output[counter], input1->get_full(), input1->get_full_element_count(), output, local_join_output[counter], local_join_output_size[counter], 1, 0, &join_tuples, iteration);
                        total_join_tuples = total_join_tuples + join_tuples;
                    }
                    else if (current_ra->get_join_input0_graph_type() == FULL && current_ra->get_join_input1_graph_type() == DELTA)
                    {
                        relation* input0 = current_ra->get_join_input0();
                        relation* output = current_ra->get_join_output();

                        join_completed = join_completed & current_ra->local_join(threshold, clique_buf_output_size[counter], &(offset[counter]), 2, clique_buf_output[counter], input0->get_full(), input0->get_full_element_count(), output, local_join_output[counter], local_join_output_size[counter], 0, 1, &join_tuples, iteration);
                        total_join_tuples = total_join_tuples + join_tuples;
                    }
            }
            counter++;
        }

        int global_synchronizer = 0;
        int synchronizer = 0;
        if (join_completed == true)
            synchronizer = 1;


        MPI_Allreduce(&synchronizer, &global_synchronizer, 1, MPI_INT, MPI_BAND, mcomm.get_comm());
        *local_join_count = total_join_tuples;
        if (global_synchronizer == 1)
        {
            u32 counter = 0;
            for (std::vector<parallel_RA*>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
            {
                parallel_RA* current_ra = *it;
                if (current_ra->get_RA_type() == JOIN)
                    delete[] clique_buf_output[counter];

                offset[counter] = 0;
                counter++;
            }

            delete[] clique_buf_output_size;
            delete[] clique_buf_output;
            return true;
        }
        else
            return false;

        //return join_completed;
        //return total_join_tuples;
    }



    int all_to_all()
    {
        u64 total_all_to_all = 0;

        u32 RA_count = RA_list.size();

        //int rank = mcomm.get_rank();
        //std::cout << "ALL TO ALL RANK " << rank << std::endl;
        u32 counter = 0;
        for (std::vector<parallel_RA*>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
        {
            parallel_RA* current_ra = *it;

            if (current_ra->get_RA_type() == COPY)
            {
                relation* output = current_ra->get_copy_output();
                current_ra->all_to_all(local_join_output[counter], local_join_output_size[counter], output);
            }
            if (current_ra->get_RA_type() == JOIN)
            {
                relation* output = current_ra->get_join_output();
                current_ra->all_to_all(local_join_output[counter], local_join_output_size[counter], output);
            }

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



    u32 local_insert_in_full(bool local_join_status)
    {
        u32 local_insert = 0;
        if (local_join_status == false)
            return local_insert;

        /*
        int counter = 0;
        for (std::vector<relation*>::iterator it = relation_list.begin() ; it != relation_list.end(); ++it)
        {
            relation* rel = *it;
            local_insert = local_insert + rel->insert_delta_in_full();
            counter++;
        }
        */

        for (std::vector<parallel_RA*>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
        {
            parallel_RA* current_ra = *it;
            if (current_ra->get_RA_type() == COPY)
            {
                relation* input = current_ra->get_copy_input();
                relation* output = current_ra->get_copy_output();

                local_insert = local_insert + input->insert_delta_in_full();
                local_insert = local_insert + output->insert_delta_in_full();

            }
            if (current_ra->get_RA_type() == JOIN)
            {
                relation* input0 = current_ra->get_join_input0();
                relation* input1 = current_ra->get_join_input1();
                //relation* output = current_ra->get_copy_output();

                local_insert = local_insert + input0->insert_delta_in_full();
                local_insert = local_insert + input1->insert_delta_in_full();
                break;
                //local_insert = local_insert + output->insert_delta_in_full();
            }
        }

        return local_insert;
    }


#if 1
    void local_insert_in_delta(bool local_join_status)
    {
        if (local_join_status == false)
            return;

        /*
        for (std::vector<relation*>::iterator it = relation_list.begin() ; it != relation_list.end(); ++it)
        {
            relation* rel = *it;
            rel->local_insert_in_delta();
        }
        */

        for (std::vector<parallel_RA*>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
        {
            parallel_RA* current_ra = *it;
            if (current_ra->get_RA_type() == COPY)
            {
                relation* input = current_ra->get_copy_input();
                relation* output = current_ra->get_copy_output();

                input->local_insert_in_delta();
                output->local_insert_in_delta();
                break;
            }
            if (current_ra->get_RA_type() == JOIN)
            {
                relation* input0 = current_ra->get_join_input0();
                relation* input1 = current_ra->get_join_input1();
                //relation* output = current_ra->get_copy_output();

                input0->local_insert_in_delta();
                input1->local_insert_in_delta();
                break;
                //output->local_insert_in_delta();
            }
        }

        return;
    }
#endif


    void execute()
    {
        //u32 iteration = 0;
        u32 clique_comm_count = 0;
        u32 local_join_count = 0;
        u32 all_to_all_count = 0;
        u32 insert_in_full_count = 0;
        //u32 insert_in_delta = 0;

        u32 running_clique_comm_count = 0;
        u32 running_local_join_count = 0;
        u32 running_all_to_all_count = 0;
        u32 running_insert_in_full_count = 0;
        //u32 running_insert_in_delta = 0;

        double clique_start = 0, clique_end = 0;
        double local_join_start = 0, local_join_end = 0;
        double all_to_all_start = 0, all_to_all_end = 0;
        double insert_full_start = 0, insert_full_end = 0;
        double insert_delta_start = 0, insert_delta_end = 0;
        double verify_start = 0, verify_end = 0;

        double iteration_time = 0;
        double running_time = 0;

        u32 RA_count = RA_list.size();
        //clique_output_queue = new std::queue<u64>*[RA_count];
        //for (u32 i = 0; i < RA_count; i++)
        //    clique_output_queue[i] = new std::queue<u64>;

        int *offset = new int[RA_count];
        for (u32 i =0; i < RA_count; i++)
            offset[i] = 0;
        bool local_join_status = true;
        int xc = 400;
        int rank = mcomm.get_rank();
        //while (xc != 0)
        int iteration = 1;
        while (true)
        {
            clique_start = MPI_Wtime();
            clique_comm_count = clique_comm(local_join_status);
            running_clique_comm_count = running_clique_comm_count + clique_comm_count;
            clique_end = MPI_Wtime();

            local_join_start = MPI_Wtime();
            local_join_status = local_join(&local_join_count, offset, iteration);
            running_local_join_count = running_local_join_count + local_join_count;
            local_join_end = MPI_Wtime();

            all_to_all_start = MPI_Wtime();
            all_to_all_count = all_to_all();
            running_all_to_all_count = running_all_to_all_count + all_to_all_count;
            all_to_all_end = MPI_Wtime();

            MPI_Barrier(mcomm.get_comm());

            insert_full_start = MPI_Wtime();
            insert_in_full_count = local_insert_in_full(local_join_status);
            running_insert_in_full_count = running_insert_in_full_count + insert_in_full_count;
            insert_full_end = MPI_Wtime();

            insert_delta_start = MPI_Wtime();
            local_insert_in_delta(local_join_status);
            //insert_in_delta = local_insert_in_delta(local_join_status);
            //running_insert_in_delta = running_insert_in_delta + insert_in_delta;
            insert_delta_end = MPI_Wtime();


            if (rank == 0)
                std::cout << "ITERATION [" << iteration <<"] " << clique_comm_count << " " << local_join_count << " " << all_to_all_count << " " << insert_in_full_count << std::endl;
            //std::cout << "RUNNING ITERATION [" << iteration <<"] " << running_clique_comm_count << " " << running_local_join_count << " " << running_all_to_all_count << " " << running_insert_in_full_count << " " << running_insert_in_delta << std::endl;


            verify_start = MPI_Wtime();
            if (check_for_fixed_point(local_join_status) == true)
            {
                verify_end = MPI_Wtime();
                iteration_time = (verify_end - verify_start) + (insert_delta_end - insert_delta_start) + (insert_full_end - insert_full_start) + (all_to_all_end - all_to_all_start) + (local_join_end - local_join_start) + (clique_end - clique_start);
                running_time = running_time + iteration_time;

                if (rank == 0)
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

            if (rank == 0)
                std::cout << " F [" << iteration << "] " << running_time << " " << iteration_time
                      << " Clique " <<  (clique_end - clique_start)
                      << " Local Join " <<  (local_join_end - local_join_start)
                      << " All to All " <<  (all_to_all_end - all_to_all_start)
                      << " Insert Full " <<  (insert_full_end - insert_full_start)
                      << " Insert delta " <<  (insert_delta_end - insert_delta_start)
                      << " Verify " <<  (verify_end - verify_start)
                      << std::endl
                      << std::endl
                      << std::endl;

            xc--;
            iteration++;
        }
        delete[] offset;


        //for (u32 i = 0; i < RA_count; i++)
        //    delete clique_output_queue[i];
        //delete[] clique_output_queue;

        print_full();
    }

};

#endif
