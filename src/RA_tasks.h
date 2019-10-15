#ifndef RAM_TASKS_H
#define RAM_TASKS_H

#include <mpi.h>
#include <vector>
#include <unordered_set>
#include "balanced_hash_relation.h"
#include "balanced_parallel_io.h"
#include "parallel_join.h"
//#include "google_btree_relation.h"
#include "btree/btree_map.h"
#include <queue>

u32 running_total = 0;
bool threshold_enabled = false;


class RAM
{
private:
    bool comm_compaction;
    u32 threshold;

    double refinement_factor;

    u32 refinement_ts;

    std::vector<parallel_RA*> RA_list;

    u64 *clique_buf_output_size;
    u64 **clique_buf_output;

    vector_buffer** local_join_output;
    int **local_join_output_size;
    int *local_cumulative_join_output_size;

    u64 *all_to_all_buffer_size;
    u64 **all_to_all_buffer;

    u64 *cumulative_all_to_all_buffer;
    int* cumulative_all_to_all_recv_process_size_array;

    mpi_comm mcomm;


public:

    RAM(mpi_comm mcomm)
    {
        this->comm_compaction = false;
        this->threshold = 1000000;
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
                    std::cout << "T Size 1 " << s0 << std::endl;
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
    }



    bool cumulative_local_join(u32* local_join_count, int* offset, int iteration)
    {
        u32 join_tuples = 0;
        u32 total_join_tuples = 0;
        u32 RA_count = RA_list.size();
        u32 nprocs = (u32)mcomm.get_nprocs();

        bool join_completed = true;

        local_join_output = new vector_buffer*[RA_count];
        local_join_output_size = new int*[RA_count];
        local_cumulative_join_output_size = new int[nprocs];
        memset(local_cumulative_join_output_size, 0, nprocs * sizeof(int));
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
                current_ra->cumulative_local_copy(input, output, local_join_output, local_join_output_size, local_cumulative_join_output_size, counter);
            }
            else if (current_ra->get_RA_type() == JOIN)
            {
                if (current_ra->get_join_input0_graph_type() == DELTA && current_ra->get_join_input1_graph_type() == DELTA)
                {
                    relation* input0 = current_ra->get_join_input0();
                    relation* output = current_ra->get_join_output();

                    join_completed = join_completed & current_ra->cumulative_local_join(threshold, clique_buf_output_size[counter], &(offset[counter]), 2, clique_buf_output[counter], input0->get_delta(), input0->get_delta_element_count(), output, local_join_output, local_join_output_size, local_cumulative_join_output_size, 0, 1, &join_tuples, counter);
                    total_join_tuples = total_join_tuples + join_tuples;
                }
                else if (current_ra->get_join_input0_graph_type() == DELTA && current_ra->get_join_input1_graph_type() == FULL)
                {
                    relation* input1 = current_ra->get_join_input1();
                    relation* output = current_ra->get_join_output();

                    join_completed = join_completed & current_ra->cumulative_local_join(threshold, clique_buf_output_size[counter], &(offset[counter]), 2, clique_buf_output[counter], input1->get_full(), input1->get_full_element_count(), output, local_join_output, local_join_output_size, local_cumulative_join_output_size, 1, 0, &join_tuples, counter);
                    total_join_tuples = total_join_tuples + join_tuples;
                }
                else if (current_ra->get_join_input0_graph_type() == FULL && current_ra->get_join_input1_graph_type() == DELTA)
                {
                    relation* input0 = current_ra->get_join_input0();
                    relation* output = current_ra->get_join_output();

                    join_completed = join_completed & current_ra->cumulative_local_join(threshold, clique_buf_output_size[counter], &(offset[counter]), 2, clique_buf_output[counter], input0->get_full(), input0->get_full_element_count(), output, local_join_output, local_join_output_size, local_cumulative_join_output_size, 0, 1, &join_tuples, counter);
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
    }



    int all_to_all()
    {
        u64 total_all_to_all = 0;

        u32 RA_count = RA_list.size();

        all_to_all_buffer_size = new u64[RA_count];
        all_to_all_buffer = new u64*[RA_count];

        u32 counter = 0;
        for (std::vector<parallel_RA*>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
        {
            parallel_RA* current_ra = *it;
            current_ra->all_to_all(local_join_output[counter], local_join_output_size[counter], &(all_to_all_buffer_size[counter]), &(all_to_all_buffer[counter]));

            delete[] local_join_output[counter];
            delete[] local_join_output_size[counter];

            counter++;
        }

        delete[] local_join_output;
        delete[] local_join_output_size;


        return total_all_to_all;
    }


    void cumulative_all_to_all(int RA_count, vector_buffer** local_join_output, int* cumulative_process_size, int** process_size, u64 **outer_hash_data, int* recv_process_size_array)
    {
        int nprocs = mcomm.get_nprocs();
        //int rank = mcomm.get_rank();

        int psize = 0;
        int psize1 = 0;
        int counter = 0;
        int *flat_process_size = new int[nprocs * RA_count];
        for (int i = 0; i < nprocs; i++)
        {
            int lc = 0;
            for (int j = 0; j < RA_count; j++)
            {
                flat_process_size[counter++] = process_size[j][i];
                psize = psize + process_size[j][i];
                psize1 = psize1 + (&local_join_output[j][i])->size;
                lc = lc + process_size[j][i];
            }
            //std::cout << "J " << j << " " << lc << std::endl;
        }

        /* This step prepares for actual data transfer */
        /* Every process sends to every other process the amount of data it is going to send */
        //int recv_process_size_array[nprocs * RA_count];
        memset(recv_process_size_array, 0, RA_count * nprocs * sizeof(int));
        MPI_Alltoall(flat_process_size, RA_count, MPI_INT, recv_process_size_array, RA_count, MPI_INT, mcomm.get_comm());

        delete[] flat_process_size;


        // PREPARING THE SEND PREFIX
        int cumulative_prefix_sum_process_size[nprocs];
        memset(cumulative_prefix_sum_process_size, 0, nprocs * sizeof(int));

        for (int i = 1; i < nprocs; i++)
            cumulative_prefix_sum_process_size[i] = cumulative_prefix_sum_process_size[i - 1] + cumulative_process_size[i - 1];

        int process_data_buffer_size = cumulative_prefix_sum_process_size[nprocs - 1] + cumulative_process_size[nprocs - 1];

        //int rank = mcomm.get_rank();
        //std::cout << rank << " XXXXXXXXXX " << psize << "   " << process_data_buffer_size << "   " << psize1 << std::endl;
        assert(psize == process_data_buffer_size);

        // PREPARING THE SEND DATA
        u64* process_data = 0;
        process_data = new u64[process_data_buffer_size];
        memset(process_data, 0, process_data_buffer_size * sizeof(u64));

        u32 boffset = 0;
        for(int i = 0; i < nprocs; i++)
        {
            for (int r = 0; r < RA_count; r++)
            {
                memcpy(process_data + boffset, (&local_join_output[r][i])->buffer, (&local_join_output[r][i])->size);
                boffset = boffset + ((&local_join_output[r][i])->size)/8;
                vector_buffer_free(&local_join_output[r][i]);
            }
        }

        // PREPARING THE RECEIVE PREFIX
        int cumulative_recv_process_size_array[nprocs];
        memset(cumulative_recv_process_size_array, 0, nprocs * sizeof(int));

        for (int k = 0; k < RA_count * nprocs; k++)
            cumulative_recv_process_size_array[k/RA_count] = cumulative_recv_process_size_array[k/RA_count] + recv_process_size_array[k];

        int cumulative_prefix_sum_recv_process_size_buffer[nprocs];
        memset(cumulative_prefix_sum_recv_process_size_buffer, 0, nprocs * sizeof(int));

        /* Sending data to all processes: What is the buffer size to allocate */
        int outer_hash_buffer_size = cumulative_recv_process_size_array[0];
        for(int i = 1; i < nprocs; i++)
        {
            cumulative_prefix_sum_recv_process_size_buffer[i] = cumulative_prefix_sum_recv_process_size_buffer[i - 1] + cumulative_recv_process_size_array[i - 1];
            outer_hash_buffer_size = outer_hash_buffer_size + cumulative_recv_process_size_array[i];
        }

        // ALLOCATE THE BUFFER
        *outer_hash_data = new u64[outer_hash_buffer_size];
        memset(*outer_hash_data, 0, outer_hash_buffer_size * sizeof(u64));


        /*
        if (rank == 0)
        {
            for (int i =0;i < nprocs; i++)
                std::cout << "AAAAAAAAAAA " <<  i << " " << cumulative_process_size[i] << " " << cumulative_prefix_sum_process_size[i] << " XX " << cumulative_recv_process_size_array[i] << " " << cumulative_prefix_sum_recv_process_size_buffer[i] << std::endl;
        }

        if (rank == 1)
        {
            for (int i =0;i < nprocs; i++)
                std::cout << "BBBBBBBBBBB " <<  i << " " << cumulative_process_size[i] << " " << cumulative_prefix_sum_process_size[i] << " XX " << cumulative_recv_process_size_array[i] << " " << cumulative_prefix_sum_recv_process_size_buffer[i] << std::endl;
        }
        */

#if 1
        // ALL TO ALL COMM
        MPI_Alltoallv(process_data, cumulative_process_size, cumulative_prefix_sum_process_size, MPI_UNSIGNED_LONG_LONG, *outer_hash_data, cumulative_recv_process_size_array, cumulative_prefix_sum_recv_process_size_buffer, MPI_UNSIGNED_LONG_LONG, mcomm.get_comm());

        delete[] process_data;
#endif
        return;
    }


    int cumulative_all_to_all()
    {
        u64 total_all_to_all = 0;
        u32 RA_count = RA_list.size();
        int nprocs = mcomm.get_nprocs();

        cumulative_all_to_all_recv_process_size_array = new int[RA_count * nprocs];


        u32 counter = 0;
        cumulative_all_to_all(RA_count, local_join_output, local_cumulative_join_output_size, local_join_output_size, &(cumulative_all_to_all_buffer), cumulative_all_to_all_recv_process_size_array);

        for (std::vector<parallel_RA*>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
        {
            //    parallel_RA* current_ra = *it;
            //    current_ra->all_to_all(local_join_output[counter], local_join_output_size[counter], &(all_to_all_buffer_size[counter]), &(all_to_all_buffer[counter]));

            delete[] local_join_output[counter];
            delete[] local_join_output_size[counter];

            counter++;
        }

        delete[] local_join_output;
        delete[] local_join_output_size;
        delete[] local_cumulative_join_output_size;



        return total_all_to_all;
    }

    int insert_in_newt1(u64 outer_hash_buffer_size, u64 *outer_hash_data, relation* output)
    {
        u32 arity = output->get_arity();
        u32 successful_insert = 0;
        u64 t[2];

        for (u64 i = 0; i < outer_hash_buffer_size; i = i + arity)
        {
            t[0] = outer_hash_data[i];
            t[1] = outer_hash_data[i + 1];

            if (output->find_in_full(t) == false)
            {
                output->insert_in_newt(t);
                successful_insert++;
            }
        }

        //int rank = mcomm.get_rank();
        //if (rank == 0)
        //{
        //int new_count = output->get_new_element_count();
        //std::cout << "Local Inserts in new " << successful_insert << " (" << new_count << ")" << std::endl;
        //}


        return successful_insert;
    }

    int cumulative_insert_in_newt( u64 *outer_hash_data, int* recv_process_size_array, std::vector<parallel_RA*> RA_list)
    {
        u32 successful_insert = 0;
        int nprocs = mcomm.get_nprocs();
        int offset = 0;
        int RA_count = RA_list.size();
        //int rank = mcomm.get_rank();

        for (int k = 0; k < RA_count * nprocs; k=k+RA_count)
        {
            successful_insert = successful_insert + insert_in_newt1(recv_process_size_array[k], (outer_hash_data) + offset, RA_list[0]->get_join_output());
            offset = offset + recv_process_size_array[k];

            successful_insert = successful_insert + insert_in_newt1(recv_process_size_array[k + 1], (outer_hash_data) + offset, RA_list[1]->get_join_output());
            offset = offset + recv_process_size_array[k + 1];

            successful_insert = successful_insert + insert_in_newt1(recv_process_size_array[k + 2], (outer_hash_data) + offset, RA_list[2]->get_join_output());
            offset = offset + recv_process_size_array[k + 2];

            successful_insert = successful_insert + insert_in_newt1(recv_process_size_array[k + 3], (outer_hash_data) + offset, RA_list[3]->get_copy_output());
            offset = offset + recv_process_size_array[k + 3];
        }

#if 0
        for (int k = 0; k < RA_count * nprocs; k++)
        {
            parallel_RA* output = RA_list[k/RA_count];
            if (output->get_RA_type() == COPY)
            {
                if (rank == 0)
                    std::cout << "C " << k << " " << k/nprocs << " " << std::endl;
                relation* copyR = output->get_copy_output();
                successful_insert = successful_insert + insert_in_newt1(recv_process_size_array[k], (outer_hash_data) + offset, copyR);
                offset = offset + recv_process_size_array[k];
            }
            else if (output->get_RA_type() == JOIN)
            {
                if (rank == 0)
                    std::cout << "J " << k << " " << k/nprocs << " " << std::endl;
                relation* copyJ = output->get_join_output();
                successful_insert = successful_insert + insert_in_newt1(recv_process_size_array[k], (outer_hash_data) + offset, copyJ);
                offset = offset + recv_process_size_array[k];
            }
        }
#endif

        return successful_insert;
    }


    void cumulative_local_insert_in_newt()
    {
        cumulative_insert_in_newt(cumulative_all_to_all_buffer, cumulative_all_to_all_recv_process_size_array, RA_list);
        delete[] cumulative_all_to_all_recv_process_size_array;
        delete[] cumulative_all_to_all_buffer;
    }


    void local_insert_in_newt()
    {
        u32 counter = 0;
        for (std::vector<parallel_RA*>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
        {
            parallel_RA* current_ra = *it;

            if (current_ra->get_RA_type() == COPY)
            {
                relation* output = current_ra->get_copy_output();
                current_ra->insert_in_newt((all_to_all_buffer_size[counter]), (all_to_all_buffer[counter]), output);
            }
            if (current_ra->get_RA_type() == JOIN)
            {
                relation* output = current_ra->get_join_output();
                current_ra->insert_in_newt((all_to_all_buffer_size[counter]), (all_to_all_buffer[counter]), output);
            }
            delete[] all_to_all_buffer[counter];
            counter++;
        }

        delete[] all_to_all_buffer_size;
        delete[] all_to_all_buffer;
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
                input->local_insert_in_delta();

                local_insert = local_insert + output->insert_delta_in_full();
                output->local_insert_in_delta();

            }
            if (current_ra->get_RA_type() == JOIN)
            {
                relation* input0 = current_ra->get_join_input0();
                relation* input1 = current_ra->get_join_input1();

                local_insert = local_insert + input0->insert_delta_in_full();
                input0->local_insert_in_delta();

                local_insert = local_insert + input1->insert_delta_in_full();
                input1->local_insert_in_delta();

                //load_balance_split_full_and_delta(input1);
                break;
            }
        }

        return local_insert;
    }


    void load_balance(bool local_join_status, float r_factor)
    {
        if (local_join_status == false)
            return;

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
            if (current_ra->get_RA_type() == JOIN)
            {
                relation* input1 = current_ra->get_join_input1();
                if (load_balance_merge_full_and_delta(input1, r_factor) == false)
                    load_balance_split_full_and_delta(input1, r_factor);
                break;
            }
        }

        return;
    }


    void initial_load_balance(float r_factor)
    {
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
            if (current_ra->get_RA_type() == JOIN)
            {
                relation* input0 = current_ra->get_join_input0();
                if (load_balance_merge_full(input0, r_factor) == false)
                    load_balance_split_full(input0, r_factor);
                break;
            }
        }

        return;
    }


    void all_to_all_comm(vector_buffer* local_join_output, int* process_size, u64 *outer_hash_buffer_size, u64 **outer_hash_data)
    {
        int nprocs = mcomm.get_nprocs();
        //int rank = mcomm.get_rank();

        /* This step prepares for actual data transfer */
        /* Every process sends to every other process the amount of data it is going to send */
        int recv_process_size_buffer[nprocs];
        memset(recv_process_size_buffer, 0, nprocs * sizeof(int));
        MPI_Alltoall(process_size, 1, MPI_INT, recv_process_size_buffer, 1, MPI_INT, mcomm.get_comm());

        int prefix_sum_process_size[nprocs];
        memset(prefix_sum_process_size, 0, nprocs * sizeof(int));

        for (int i = 1; i < nprocs; i++)
            prefix_sum_process_size[i] = prefix_sum_process_size[i - 1] + process_size[i - 1];

        int process_data_buffer_size = prefix_sum_process_size[nprocs - 1] + process_size[nprocs - 1];

        u64* process_data = 0;
        process_data = new u64[process_data_buffer_size];
        memset(process_data, 0, process_data_buffer_size * sizeof(u64));

        for(int i = 0; i < nprocs; i++)
        {
            memcpy(process_data + prefix_sum_process_size[i], (&local_join_output[i])->buffer, (&local_join_output[i])->size);
            vector_buffer_free(&local_join_output[i]);
        }

        int prefix_sum_recv_process_size_buffer[nprocs];
        memset(prefix_sum_recv_process_size_buffer, 0, nprocs * sizeof(int));

        /* Sending data to all processes: What is the buffer size to allocate */
        *outer_hash_buffer_size = recv_process_size_buffer[0];
        for(int i = 1; i < nprocs; i++)
        {
            prefix_sum_recv_process_size_buffer[i] = prefix_sum_recv_process_size_buffer[i - 1] + recv_process_size_buffer[i - 1];
            *outer_hash_buffer_size = *outer_hash_buffer_size + recv_process_size_buffer[i];
        }

        *outer_hash_data = new u64[*outer_hash_buffer_size];
        memset(*outer_hash_data, 0, *outer_hash_buffer_size * sizeof(u64));


        MPI_Alltoallv(process_data, process_size, prefix_sum_process_size, MPI_UNSIGNED_LONG_LONG, *outer_hash_data, recv_process_size_buffer, prefix_sum_recv_process_size_buffer, MPI_UNSIGNED_LONG_LONG, mcomm.get_comm());

        delete[] process_data;

        return;
    }

    bool load_balance_split_full_and_delta(relation* rel, float rf)
    {
        u32 rank = (u32)mcomm.get_rank();
        int nprocs = mcomm.get_nprocs();
        int buckets = mcomm.get_number_of_buckets();
        u32* bucket_map = rel->get_bucket_map();
        u32* sub_bucket_count = rel->get_sub_bucket_count();
        google_relation* full = rel->get_full();
        google_relation* delta = rel->get_delta();
        u32** full_sub_bucket_size = rel->get_full_sub_bucket_element_count();
        u32** delta_sub_bucket_size = rel->get_delta_sub_bucket_element_count();
        u32** newt_sub_bucket_size = rel->get_new_sub_bucket_element_count();
        u32** sub_bucket_rank = rel->get_sub_bucket_rank();
        int** distinct_sub_bucket_rank = rel->get_distinct_sub_bucket_rank();
        int* distinct_sub_bucket_rank_count = rel->get_distinct_sub_bucket_rank_count();
        int arity = rel->get_arity();

        int min_sub_bucket_size = INT_MAX;
        int *max_sub_bucket_size = new int[buckets];
        memset(max_sub_bucket_size, 0, buckets * sizeof(int));


#if 0
        u64 T_tuple_count_before = 0;
        u64 T_tuple_count_before_global = 0;

        for (int bk = 0; bk < buckets; bk++)
        {
            for (u32 j = 0; j < sub_bucket_count[bk]; j++)
                T_tuple_count_before = T_tuple_count_before + full_sub_bucket_size[bk][j];
        }
        MPI_Allreduce(&T_tuple_count_before, &T_tuple_count_before_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        u64 dT_tuple_count_before = 0;
        u64 dT_tuple_count_before_global = 0;

        for (int bk = 0; bk < buckets; bk++)
        {
            for (u32 j = 0; j < sub_bucket_count[bk]; j++)
                dT_tuple_count_before = dT_tuple_count_before + delta_sub_bucket_size[bk][j];
        }
        MPI_Allreduce(&dT_tuple_count_before, &dT_tuple_count_before_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        //if (rank == 0)
        //    std::cout << "[T] Before Load balancing " << T_tuple_count_before_global << std::endl;
#endif

        u32 global_total_sub_bucket_size = 0;
        u32 total_sub_bucket_size = 0;
        u32 total_sub_bucket_count = 0;
        for (int i = 0; i < buckets; i++)
        {
            total_sub_bucket_count = total_sub_bucket_count + sub_bucket_count[i];
            if (bucket_map[i] == 1)
            {
                for (u32 j = 0; j < sub_bucket_count[i]; j++)
                {
                    if (full_sub_bucket_size[i][j] != 0)
                    {
                        if ((int)full_sub_bucket_size[i][j] > max_sub_bucket_size[i])
                            max_sub_bucket_size[i] = full_sub_bucket_size[i][j];

                        if ((int)full_sub_bucket_size[i][j] < min_sub_bucket_size)
                            min_sub_bucket_size = full_sub_bucket_size[i][j];

                        total_sub_bucket_size = total_sub_bucket_size + full_sub_bucket_size[i][j];
                    }
                }
            }
        }

        int *global_max = new int[buckets];
        memset(global_max, 0, buckets * sizeof(int));

        int global_min;
        int average_sub_bucket_size;
        MPI_Allreduce(max_sub_bucket_size, global_max, buckets, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&min_sub_bucket_size, &global_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&total_sub_bucket_size, &global_total_sub_bucket_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        delete[] max_sub_bucket_size;

        average_sub_bucket_size = global_total_sub_bucket_size / total_sub_bucket_count;

        u32 global_new_sub_bucket[buckets];
        memcpy(global_new_sub_bucket, sub_bucket_count, buckets * sizeof(u32));

        int count = 0;
        int old_total_sub_buckets = 0;
        int new_total_sub_buckets = 0;
        int global_global_max = 0;
        int average_global_max = 0;
        for (int i = 0; i < buckets; i++)
        {
            if (global_new_sub_bucket[i] >= (u32)nprocs)
                continue;
            if (global_max[i] > average_sub_bucket_size * rf * 0.8)
                //if (global_max[i] > global_min * rf * 0.8)
                //if (global_max[i] > global_min * rf)
            {
                //if (rank == 0)
                //    std::cout << "XXXX " << global_max[i] << " " << average_sub_bucket_size * rf * 0.8 << std::endl;
                global_new_sub_bucket[i] = global_new_sub_bucket[i] * rf;
                count++;
            }

            average_global_max = average_global_max + global_max[i];

            if (global_global_max < global_max[i])
                global_global_max = global_max[i];


            old_total_sub_buckets = old_total_sub_buckets + sub_bucket_count[i];
            new_total_sub_buckets = new_total_sub_buckets + global_new_sub_bucket[i];
        }
        delete[] global_max;

        if (count == 0)
        {
            if (rank == 0)
                std::cout << "[NO] RF " << rf << " Bucket Split -- Global Min " << global_min << " Average bucket size " << average_sub_bucket_size << " Global Max [" << global_global_max << " " << average_global_max/buckets  << "] OLD " << old_total_sub_buckets << " NEW " << new_total_sub_buckets << std::endl;
            return false;
        }
        else if (count != 0)
        {
            if (rank == 0)
                std::cout << "[YES] RF " << rf << " Bucket Split -- Global Min " << global_min << " Average bucket size " << average_sub_bucket_size << " Global Max [" << global_global_max << " " << average_global_max/buckets  << "] OLD " << old_total_sub_buckets << " NEW " << new_total_sub_buckets << std::endl;
        }


        int rcount = sub_bucket_rank[buckets-1][sub_bucket_count[buckets-1] - 1] + 1;

        for (int b = 0; b < buckets; b++)
        {
            if (sub_bucket_count[b] < global_new_sub_bucket[b])
            {
                bucket_map[b] = 0;

                std::unordered_set<int> distinct_t_ranks;
                u32* temp = new u32[sub_bucket_count[b]];
                memcpy(temp, sub_bucket_rank[b], sizeof(int) * sub_bucket_count[b]);

                delete[] sub_bucket_rank[b];
                sub_bucket_rank[b] = new u32[global_new_sub_bucket[b]];

                memcpy(sub_bucket_rank[b], temp, sizeof(int) * sub_bucket_count[b]);
                delete[] temp;

                for (u64 x = 0; x < global_new_sub_bucket[b]; x++)
                {
                    if (x >= sub_bucket_count[b])
                    {
                        sub_bucket_rank[b][x] = rcount % nprocs;
                        rcount++;
                    }

                    if (sub_bucket_rank[b][x] == rank)
                        bucket_map[b] = 1;

                    distinct_t_ranks.insert(sub_bucket_rank[b][x]);

                }

                delete[] distinct_sub_bucket_rank[b];
                distinct_sub_bucket_rank[b] = new int[distinct_t_ranks.size()];
                u32 x = 0;
                for ( auto it = distinct_t_ranks.begin(); it != distinct_t_ranks.end(); ++it )
                {
                    distinct_sub_bucket_rank[b][x] = *it;
                    x++;
                }
                distinct_sub_bucket_rank_count[b] = x;
            }

            //if (rank == 0)
            //    std::cout << "GNSB " << b << " " << global_new_sub_bucket[b] << std::endl;
        }


        u32* old_sub_bucket_count = new u32[buckets];
        memcpy(old_sub_bucket_count, sub_bucket_count, buckets * sizeof(u32));
        memcpy(sub_bucket_count, global_new_sub_bucket, sizeof(u32) * buckets);


        int process_size[nprocs];
        memset(process_size, 0, nprocs * sizeof(int));

        vector_buffer* process_data_vector = (vector_buffer*)malloc(sizeof(vector_buffer) * nprocs);
        for (int i = 0; i < nprocs; ++i) {
            process_data_vector[i] = vector_buffer_create_empty();
        }

        int prefix_sum_process_size[nprocs];
        memset(prefix_sum_process_size, 0, nprocs * sizeof(int));

        int recv_process_size_buffer[nprocs];
        memset(recv_process_size_buffer, 0, nprocs * sizeof(int));

        int process_size_dt[nprocs];
        memset(process_size_dt, 0, nprocs * sizeof(int));

        vector_buffer* process_data_vector_dt1 = (vector_buffer*)malloc(sizeof(vector_buffer) * nprocs);
        for (int i = 0; i < nprocs; ++i) {
            process_data_vector_dt1[i] = vector_buffer_create_empty();
        }

        int prefix_sum_process_size_dt[nprocs];
        memset(prefix_sum_process_size_dt, 0, nprocs * sizeof(int));

        int recv_process_size_buffer_dt[nprocs];
        memset(recv_process_size_buffer_dt, 0, nprocs * sizeof(int));

        u64 val[2] = {0, 0};
        for (int bk = 0; bk < buckets; bk++)
        {
            if (old_sub_bucket_count[bk] != global_new_sub_bucket[bk])
            {
                delete[] full_sub_bucket_size[bk];
                full_sub_bucket_size[bk] = new u32[global_new_sub_bucket[bk]];
                memset(full_sub_bucket_size[bk], 0, sizeof(u32) * global_new_sub_bucket[bk]);

                delete[] delta_sub_bucket_size[bk];
                delta_sub_bucket_size[bk] = new u32[global_new_sub_bucket[bk]];
                memset(delta_sub_bucket_size[bk], 0, sizeof(u32) * global_new_sub_bucket[bk]);

                delete[] newt_sub_bucket_size[bk];
                newt_sub_bucket_size[bk] = new u32[global_new_sub_bucket[bk]];
                memset(newt_sub_bucket_size[bk], 0, sizeof(u32) * global_new_sub_bucket[bk]);

                for (auto it = full[bk].begin(); it != full[bk].end(); it++)
                {
                    Map0* it2 = it->second;
                    for (auto dit2 = it2->begin(); dit2 != it2->end(); dit2++)
                    {
                        val[0] = it->first;
                        val[1] = dit2->first;
                        uint64_t bucket_id = hash_function(val[0]) % buckets;
                        uint64_t sub_bucket_id = hash_function(val[1]) % global_new_sub_bucket[bucket_id];

                        int index = sub_bucket_rank[bucket_id][sub_bucket_id];
                        process_size[index] = process_size[index] + arity;


                        vector_buffer_append(&process_data_vector[index],(unsigned char *) val, sizeof(u64) * arity);

                        //if (index == 0)
                        //std::cout << "BI and SBI " << val[0] << " " << val[1] << " " << hash_function(val[0]) << " " << bucket_id << " " << sub_bucket_id << " rank " << index << std::endl;
                    }
                }

                for(google_relation::iterator iy2 = full[bk].begin(); iy2 != full[bk].end(); iy2++)
                    delete (iy2->second);

                full[bk].clear();

                for (auto it = delta[bk].begin(); it != delta[bk].end(); it++)
                {
                    Map0* it2 = it->second;
                    for (auto dit2 = it2->begin(); dit2 != it2->end(); dit2++)
                    {
                        val[0] = it->first;
                        val[1] = dit2->first;
                        uint64_t bucket_id = hash_function(val[0]) % buckets;
                        uint64_t sub_bucket_id = hash_function(val[1]) % global_new_sub_bucket[bucket_id];

                        int index = sub_bucket_rank[bucket_id][sub_bucket_id];
                        process_size_dt[index] = process_size_dt[index] + arity;


                        vector_buffer_append(&process_data_vector_dt1[index],(unsigned char *) val, sizeof(u64)*arity);
                    }
                }

                for(google_relation::iterator iy2 = delta[bk].begin(); iy2 != delta[bk].end(); iy2++)
                    delete (iy2->second);

                delta[bk].clear();
            }
        }

        u64 outer_hash_buffer_size = 0;
        u64* outer_hash_data;
        all_to_all_comm(process_data_vector, process_size, &outer_hash_buffer_size, &outer_hash_data);
        free(process_data_vector);

        //if (rank == 0)
        //    std::cout << "outer_hash_buffer_size " << outer_hash_buffer_size <<std::endl;

        u64 t[2];
        for (u64 in = 0; in < outer_hash_buffer_size; in = in + arity)
        {
            t[0] = outer_hash_data[in];
            t[1] = outer_hash_data[in + 1];
            rel->insert_in_full(t);
        }

        delete[] outer_hash_data;


        u64 dt_outer_hash_buffer_size = 0;
        u64* dt_outer_hash_data;
        all_to_all_comm(process_data_vector_dt1, process_size_dt, &dt_outer_hash_buffer_size, &dt_outer_hash_data);
        free(process_data_vector_dt1);

        for (u64 in = 0; in < dt_outer_hash_buffer_size; in = in + 2)
        {
            t[0] = dt_outer_hash_data[in];
            t[1] = dt_outer_hash_data[in + 1];
            rel->insert_in_delta(t);
        }


        delete[] old_sub_bucket_count;
        delete[] dt_outer_hash_data;

#if 0
        u64 T_tuple_count_after = 0;
        u64 T_tuple_count_after_global = 0;
        for (int bk = 0; bk < buckets; bk++)
        {
            for (u32 j = 0; j < global_new_sub_bucket[bk]; j++)
                T_tuple_count_after = T_tuple_count_after + full_sub_bucket_size[bk][j];
        }
        MPI_Allreduce(&T_tuple_count_after, &T_tuple_count_after_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        //if (rank == 0)
        //    std::cout << "[T] After Load balancing " << T_tuple_count_before_global << " " << T_tuple_count_after_global << std::endl;
        assert(T_tuple_count_before_global == T_tuple_count_after_global);


        u64 dT_tuple_count_after = 0;
        u64 dT_tuple_count_after_global = 0;
        for (int bk = 0; bk < buckets; bk++)
        {
            for (u32 j = 0; j < global_new_sub_bucket[bk]; j++)
                dT_tuple_count_after = dT_tuple_count_after + delta_sub_bucket_size[bk][j];
        }
        MPI_Allreduce(&dT_tuple_count_after, &dT_tuple_count_after_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        //if (rank == 0)
        //    std::cout << "[dT] After Load balancing " << dT_tuple_count_before_global << " " << dT_tuple_count_after_global << std::endl;
        assert(dT_tuple_count_before_global == dT_tuple_count_after_global);
#endif

        return true;
    }


    bool load_balance_split_full(relation* rel, float rf)
    {
        u32 rank = (u32)mcomm.get_rank();
        int nprocs = mcomm.get_nprocs();
        int buckets = mcomm.get_number_of_buckets();
        u32* bucket_map = rel->get_bucket_map();
        u32* sub_bucket_count = rel->get_sub_bucket_count();
        google_relation* full = rel->get_full();
        u32** full_sub_bucket_size = rel->get_full_sub_bucket_element_count();
        u32** sub_bucket_rank = rel->get_sub_bucket_rank();
        int** distinct_sub_bucket_rank = rel->get_distinct_sub_bucket_rank();
        int* distinct_sub_bucket_rank_count = rel->get_distinct_sub_bucket_rank_count();
        int arity = rel->get_arity();

        int min_sub_bucket_size = INT_MAX;
        int *max_sub_bucket_size = new int[buckets];
        memset(max_sub_bucket_size, 0, buckets * sizeof(int));

        u32 global_total_sub_bucket_size = 0;
        u32 total_sub_bucket_size = 0;
        u32 total_sub_bucket_count = 0;
        for (int i = 0; i < buckets; i++)
        {
            total_sub_bucket_count = total_sub_bucket_count + sub_bucket_count[i];
            if (bucket_map[i] == 1)
            {
                for (u32 j = 0; j < sub_bucket_count[i]; j++)
                {
                    if (full_sub_bucket_size[i][j] != 0)
                    {
                        if ((int)full_sub_bucket_size[i][j] > max_sub_bucket_size[i])
                            max_sub_bucket_size[i] = full_sub_bucket_size[i][j];

                        if ((int)full_sub_bucket_size[i][j] < min_sub_bucket_size)
                            min_sub_bucket_size = full_sub_bucket_size[i][j];

                        total_sub_bucket_size = total_sub_bucket_size + full_sub_bucket_size[i][j];
                    }
                }
            }
        }

        int *global_max = new int[buckets];
        memset(global_max, 0, buckets * sizeof(int));

        int global_min;
        int average_sub_bucket_size;
        MPI_Allreduce(max_sub_bucket_size, global_max, buckets, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&min_sub_bucket_size, &global_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&total_sub_bucket_size, &global_total_sub_bucket_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        delete[] max_sub_bucket_size;

        average_sub_bucket_size = global_total_sub_bucket_size / total_sub_bucket_count;

        u32 global_new_sub_bucket[buckets];
        memcpy(global_new_sub_bucket, sub_bucket_count, buckets * sizeof(u32));

        int count = 0;
        int old_total_sub_buckets = 0;
        int new_total_sub_buckets = 0;
        int global_global_max = 0;
        int average_global_max = 0;
        for (int i = 0; i < buckets; i++)
        {
            if (global_new_sub_bucket[i] >= (u32)nprocs)
                continue;
            if (global_max[i] > average_sub_bucket_size * rf * 0.8)
                //if (global_max[i] > global_min * rf * 0.8)
                //if (global_max[i] > global_min * rf)
            {
                //if (rank == 0)
                //    std::cout << "XXXX " << global_max[i] << " " << average_sub_bucket_size * rf * 0.8 << std::endl;
                global_new_sub_bucket[i] = global_new_sub_bucket[i] * rf;
                count++;
            }

            average_global_max = average_global_max + global_max[i];

            if (global_global_max < global_max[i])
                global_global_max = global_max[i];


            old_total_sub_buckets = old_total_sub_buckets + sub_bucket_count[i];
            new_total_sub_buckets = new_total_sub_buckets + global_new_sub_bucket[i];
        }
        delete[] global_max;

        if (count == 0)
        {
            if (rank == 0)
                std::cout << "[NO] RF " << rf << " Bucket Split -- Global Min " << global_min << " Average bucket size " << average_sub_bucket_size << " Global Max [" << global_global_max << " " << average_global_max/buckets  << "] OLD " << old_total_sub_buckets << " NEW " << new_total_sub_buckets << std::endl;
            return false;
        }
        else if (count != 0)
        {
            if (rank == 0)
                std::cout << "[YES] RF " << rf << " Bucket Split -- Global Min " << global_min << " Average bucket size " << average_sub_bucket_size << " Global Max [" << global_global_max << " " << average_global_max/buckets  << "] OLD " << old_total_sub_buckets << " NEW " << new_total_sub_buckets << std::endl;
        }


        int rcount = sub_bucket_rank[buckets-1][sub_bucket_count[buckets-1] - 1] + 1;

        for (int b = 0; b < buckets; b++)
        {
            if (sub_bucket_count[b] < global_new_sub_bucket[b])
            {
                bucket_map[b] = 0;

                std::unordered_set<int> distinct_t_ranks;
                u32* temp = new u32[sub_bucket_count[b]];
                memcpy(temp, sub_bucket_rank[b], sizeof(int) * sub_bucket_count[b]);

                delete[] sub_bucket_rank[b];
                sub_bucket_rank[b] = new u32[global_new_sub_bucket[b]];

                memcpy(sub_bucket_rank[b], temp, sizeof(int) * sub_bucket_count[b]);
                delete[] temp;

                for (u64 x = 0; x < global_new_sub_bucket[b]; x++)
                {
                    if (x >= sub_bucket_count[b])
                    {
                        sub_bucket_rank[b][x] = rcount % nprocs;
                        rcount++;
                    }

                    if (sub_bucket_rank[b][x] == rank)
                        bucket_map[b] = 1;

                    distinct_t_ranks.insert(sub_bucket_rank[b][x]);

                }

                delete[] distinct_sub_bucket_rank[b];
                distinct_sub_bucket_rank[b] = new int[distinct_t_ranks.size()];
                u32 x = 0;
                for ( auto it = distinct_t_ranks.begin(); it != distinct_t_ranks.end(); ++it )
                {
                    distinct_sub_bucket_rank[b][x] = *it;
                    x++;
                }
                distinct_sub_bucket_rank_count[b] = x;
            }

            //if (rank == 0)
            //    std::cout << "GNSB " << b << " " << global_new_sub_bucket[b] << std::endl;
        }


        u32* old_sub_bucket_count = new u32[buckets];
        memcpy(old_sub_bucket_count, sub_bucket_count, buckets * sizeof(u32));
        memcpy(sub_bucket_count, global_new_sub_bucket, sizeof(u32) * buckets);


        int process_size[nprocs];
        memset(process_size, 0, nprocs * sizeof(int));

        vector_buffer* process_data_vector = (vector_buffer*)malloc(sizeof(vector_buffer) * nprocs);
        for (int i = 0; i < nprocs; ++i) {
            process_data_vector[i] = vector_buffer_create_empty();
        }

        int prefix_sum_process_size[nprocs];
        memset(prefix_sum_process_size, 0, nprocs * sizeof(int));

        int recv_process_size_buffer[nprocs];
        memset(recv_process_size_buffer, 0, nprocs * sizeof(int));

        int process_size_dt[nprocs];
        memset(process_size_dt, 0, nprocs * sizeof(int));

        vector_buffer* process_data_vector_dt1 = (vector_buffer*)malloc(sizeof(vector_buffer) * nprocs);
        for (int i = 0; i < nprocs; ++i) {
            process_data_vector_dt1[i] = vector_buffer_create_empty();
        }

        int prefix_sum_process_size_dt[nprocs];
        memset(prefix_sum_process_size_dt, 0, nprocs * sizeof(int));

        int recv_process_size_buffer_dt[nprocs];
        memset(recv_process_size_buffer_dt, 0, nprocs * sizeof(int));

        u64 val[2] = {0, 0};
        for (int bk = 0; bk < buckets; bk++)
        {
            if (old_sub_bucket_count[bk] != global_new_sub_bucket[bk])
            {
                delete[] full_sub_bucket_size[bk];
                full_sub_bucket_size[bk] = new u32[global_new_sub_bucket[bk]];
                memset(full_sub_bucket_size[bk], 0, sizeof(u32) * global_new_sub_bucket[bk]);

                for (auto it = full[bk].begin(); it != full[bk].end(); it++)
                {
                    Map0* it2 = it->second;
                    for (auto dit2 = it2->begin(); dit2 != it2->end(); dit2++)
                    {
                        val[0] = it->first;
                        val[1] = dit2->first;
                        uint64_t bucket_id = hash_function(val[0]) % buckets;
                        uint64_t sub_bucket_id = hash_function(val[1]) % global_new_sub_bucket[bucket_id];

                        int index = sub_bucket_rank[bucket_id][sub_bucket_id];
                        process_size[index] = process_size[index] + arity;


                        vector_buffer_append(&process_data_vector[index],(unsigned char *) val, sizeof(u64) * arity);

                        //if (index == 0)
                        //std::cout << "BI and SBI " << val[0] << " " << val[1] << " " << hash_function(val[0]) << " " << bucket_id << " " << sub_bucket_id << " rank " << index << std::endl;
                    }
                }

                for(google_relation::iterator iy2 = full[bk].begin(); iy2 != full[bk].end(); iy2++)
                    delete (iy2->second);

                full[bk].clear();
            }
        }

        u64 outer_hash_buffer_size = 0;
        u64* outer_hash_data;
        all_to_all_comm(process_data_vector, process_size, &outer_hash_buffer_size, &outer_hash_data);
        free(process_data_vector);

        u64 t[2];
        for (u64 in = 0; in < outer_hash_buffer_size; in = in + arity)
        {
            t[0] = outer_hash_data[in];
            t[1] = outer_hash_data[in + 1];
            rel->insert_in_full(t);
        }

        delete[] outer_hash_data;
        delete[] old_sub_bucket_count;


        return true;
    }


    bool load_balance_merge_full_and_delta(relation* rel, float rf)
    {
        u32 rank = (u32)mcomm.get_rank();
        int nprocs = mcomm.get_nprocs();
        int buckets = mcomm.get_number_of_buckets();
        u32* bucket_map = rel->get_bucket_map();
        u32* sub_bucket_count = rel->get_sub_bucket_count();
        google_relation* full = rel->get_full();
        google_relation* delta = rel->get_delta();
        u32** full_sub_bucket_size = rel->get_full_sub_bucket_element_count();
        u32** delta_sub_bucket_size = rel->get_delta_sub_bucket_element_count();
        u32** newt_sub_bucket_size = rel->get_new_sub_bucket_element_count();
        u32** sub_bucket_rank = rel->get_sub_bucket_rank();
        int** distinct_sub_bucket_rank = rel->get_distinct_sub_bucket_rank();
        int* distinct_sub_bucket_rank_count = rel->get_distinct_sub_bucket_rank_count();
        int arity = rel->get_arity();

        u32 global_new_sub_bucket[buckets];
        memcpy(global_new_sub_bucket, sub_bucket_count, buckets * sizeof(u32));

        int mcount = 0;
        u32 new_sub_bucket_count = 0;
        u32 total_old_sub_bucket_count = 0;
        for (int i = 0; i < buckets; i++)
        {
            total_old_sub_bucket_count = total_old_sub_bucket_count + sub_bucket_count[i];
            if (sub_bucket_count[i] >= rf)
                mcount++;
        }


        if (mcount > 0.6* buckets)
        {
            for (int i = 0; i < buckets; i++)
            {
                global_new_sub_bucket[i] = global_new_sub_bucket[i] / rf;
                if (global_new_sub_bucket[i] == 0)
                    global_new_sub_bucket[i] = 1;

                new_sub_bucket_count = new_sub_bucket_count + global_new_sub_bucket[i];
            }
            if (rank == 0)
                std::cout << "[YES] Bucket Consolidation [" << mcount << " " << 0.6 * buckets << "] Old Sub Bucket count " << total_old_sub_bucket_count << " New Sub Bucket count " << new_sub_bucket_count << std::endl;
        }
        else
        {
            if (rank == 0)
                std::cout << "[NO] Bucket Consolidation [" << mcount << " " << 0.6 * buckets << "] Old Sub Bucket count " << total_old_sub_bucket_count << std::endl;
            return false;
        }

#if 0
        u64 T_tuple_count_before = 0;
        u64 T_tuple_count_before_global = 0;

        for (int bk = 0; bk < buckets; bk++)
        {
            for (u32 j = 0; j < sub_bucket_count[bk]; j++)
                T_tuple_count_before = T_tuple_count_before + full_sub_bucket_size[bk][j];
        }
        MPI_Allreduce(&T_tuple_count_before, &T_tuple_count_before_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        u64 dT_tuple_count_before = 0;
        u64 dT_tuple_count_before_global = 0;

        for (int bk = 0; bk < buckets; bk++)
        {
            for (u32 j = 0; j < sub_bucket_count[bk]; j++)
                dT_tuple_count_before = dT_tuple_count_before + delta_sub_bucket_size[bk][j];
        }
        MPI_Allreduce(&dT_tuple_count_before, &dT_tuple_count_before_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        //if (rank == 0)
        //    std::cout << "[T] Before Load balancing " << T_tuple_count_before_global << std::endl;
#endif

        int rcount = sub_bucket_rank[buckets-1][sub_bucket_count[buckets-1] - 1] + 1;
        for (int b = 0; b < buckets; b++)
        {
            if (sub_bucket_count[b] > global_new_sub_bucket[b])
            {
                bucket_map[b] = 0;

                std::unordered_set<int> distinct_t_ranks;
                int* temp = new int[sub_bucket_count[b]];
                memcpy(temp, sub_bucket_rank[b], sizeof(int) * sub_bucket_count[b]);

                delete[] sub_bucket_rank[b];
                sub_bucket_rank[b] = new u32[global_new_sub_bucket[b]];

                memcpy(sub_bucket_rank[b], temp, sizeof(int) * global_new_sub_bucket[b]);
                delete[] temp;

                for (u64 x = 0; x < global_new_sub_bucket[b]; x++)
                {
                    if (sub_bucket_rank[b][x] == rank)
                        bucket_map[b] = 1;

                    distinct_t_ranks.insert(sub_bucket_rank[b][x]);
                    rcount++;
                }

                delete[] distinct_sub_bucket_rank[b];
                distinct_sub_bucket_rank[b] = new int[distinct_t_ranks.size()];
                u32 x = 0;
                for ( auto it = distinct_t_ranks.begin(); it != distinct_t_ranks.end(); ++it )
                {
                    distinct_sub_bucket_rank[b][x] = *it;
                    x++;
                }
                distinct_sub_bucket_rank_count[b] = x;
            }

            //if (rank == 0)
            //    std::cout << "GNSB " << b << " " << global_new_sub_bucket[b] << std::endl;
        }


        u32* old_sub_bucket_count = new u32[buckets];
        memcpy(old_sub_bucket_count, sub_bucket_count, buckets * sizeof(u32));
        memcpy(sub_bucket_count, global_new_sub_bucket, sizeof(u32) * buckets);


        int process_size[nprocs];
        memset(process_size, 0, nprocs * sizeof(int));

        vector_buffer* process_data_vector = (vector_buffer*)malloc(sizeof(vector_buffer) * nprocs);
        for (int i = 0; i < nprocs; ++i) {
            process_data_vector[i] = vector_buffer_create_empty();
        }

        int prefix_sum_process_size[nprocs];
        memset(prefix_sum_process_size, 0, nprocs * sizeof(int));

        int recv_process_size_buffer[nprocs];
        memset(recv_process_size_buffer, 0, nprocs * sizeof(int));

        int process_size_dt[nprocs];
        memset(process_size_dt, 0, nprocs * sizeof(int));

        vector_buffer* process_data_vector_dt1 = (vector_buffer*)malloc(sizeof(vector_buffer) * nprocs);
        for (int i = 0; i < nprocs; ++i) {
            process_data_vector_dt1[i] = vector_buffer_create_empty();
        }

        int prefix_sum_process_size_dt[nprocs];
        memset(prefix_sum_process_size_dt, 0, nprocs * sizeof(int));

        int recv_process_size_buffer_dt[nprocs];
        memset(recv_process_size_buffer_dt, 0, nprocs * sizeof(int));

        u64 val[2] = {0, 0};
        for (int bk = 0; bk < buckets; bk++)
        {
            if (old_sub_bucket_count[bk] != global_new_sub_bucket[bk])
            {
                delete[] full_sub_bucket_size[bk];
                full_sub_bucket_size[bk] = new u32[global_new_sub_bucket[bk]];
                memset(full_sub_bucket_size[bk], 0, sizeof(u32) * global_new_sub_bucket[bk]);

                delete[] delta_sub_bucket_size[bk];
                delta_sub_bucket_size[bk] = new u32[global_new_sub_bucket[bk]];
                memset(delta_sub_bucket_size[bk], 0, sizeof(u32) * global_new_sub_bucket[bk]);

                delete[] newt_sub_bucket_size[bk];
                newt_sub_bucket_size[bk] = new u32[global_new_sub_bucket[bk]];
                memset(newt_sub_bucket_size[bk], 0, sizeof(u32) * global_new_sub_bucket[bk]);

                for (auto it = full[bk].begin(); it != full[bk].end(); it++)
                {
                    Map0* it2 = it->second;
                    for (auto dit2 = it2->begin(); dit2 != it2->end(); dit2++)
                    {
                        val[0] = it->first;
                        val[1] = dit2->first;
                        uint64_t bucket_id = hash_function(val[0]) % buckets;
                        uint64_t sub_bucket_id = hash_function(val[1]) % global_new_sub_bucket[bucket_id];

                        int index = sub_bucket_rank[bucket_id][sub_bucket_id];
                        process_size[index] = process_size[index] + arity;


                        vector_buffer_append(&process_data_vector[index],(unsigned char *) val, sizeof(u64) * arity);

                        //if (index == 0)
                        //std::cout << "BI and SBI " << val[0] << " " << val[1] << " " << hash_function(val[0]) << " " << bucket_id << " " << sub_bucket_id << " rank " << index << std::endl;
                    }
                }

                for(google_relation::iterator iy2 = full[bk].begin(); iy2 != full[bk].end(); iy2++)
                    delete (iy2->second);

                full[bk].clear();

                for (auto it = delta[bk].begin(); it != delta[bk].end(); it++)
                {
                    Map0* it2 = it->second;
                    for (auto dit2 = it2->begin(); dit2 != it2->end(); dit2++)
                    {
                        val[0] = it->first;
                        val[1] = dit2->first;
                        uint64_t bucket_id = hash_function(val[0]) % buckets;
                        uint64_t sub_bucket_id = hash_function(val[1]) % global_new_sub_bucket[bucket_id];

                        int index = sub_bucket_rank[bucket_id][sub_bucket_id];
                        process_size_dt[index] = process_size_dt[index] + arity;


                        vector_buffer_append(&process_data_vector_dt1[index],(unsigned char *) val, sizeof(u64)*arity);
                    }
                }

                for(google_relation::iterator iy2 = delta[bk].begin(); iy2 != delta[bk].end(); iy2++)
                    delete (iy2->second);

                delta[bk].clear();
            }
        }

        u64 outer_hash_buffer_size = 0;
        u64* outer_hash_data;
        all_to_all_comm(process_data_vector, process_size, &outer_hash_buffer_size, &outer_hash_data);

        free (process_data_vector);

        //if (rank == 0)
        //    std::cout << "outer_hash_buffer_size " << outer_hash_buffer_size <<std::endl;

        u64 t[2];
        for (u64 in = 0; in < outer_hash_buffer_size; in = in + arity)
        {
            t[0] = outer_hash_data[in];
            t[1] = outer_hash_data[in + 1];
            rel->insert_in_full(t);
        }
        delete[] outer_hash_data;

        u64 dt_outer_hash_buffer_size = 0;
        u64* dt_outer_hash_data;
        all_to_all_comm(process_data_vector_dt1, process_size_dt, &dt_outer_hash_buffer_size, &dt_outer_hash_data);

        free (process_data_vector_dt1);

        for (u64 in = 0; in < dt_outer_hash_buffer_size; in = in + 2)
        {
            t[0] = dt_outer_hash_data[in];
            t[1] = dt_outer_hash_data[in + 1];
            rel->insert_in_delta(t);
        }

        delete[] old_sub_bucket_count;
        delete[] dt_outer_hash_data;

#if 0
        u64 T_tuple_count_after = 0;
        u64 T_tuple_count_after_global = 0;
        for (int bk = 0; bk < buckets; bk++)
        {
            for (u32 j = 0; j < global_new_sub_bucket[bk]; j++)
                T_tuple_count_after = T_tuple_count_after + full_sub_bucket_size[bk][j];
        }
        MPI_Allreduce(&T_tuple_count_after, &T_tuple_count_after_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        //if (rank == 0)
        //    std::cout << "[T] After Load balancing " << T_tuple_count_before_global << " " << T_tuple_count_after_global << std::endl;
        assert(T_tuple_count_before_global == T_tuple_count_after_global);


        u64 dT_tuple_count_after = 0;
        u64 dT_tuple_count_after_global = 0;
        for (int bk = 0; bk < buckets; bk++)
        {
            for (u32 j = 0; j < global_new_sub_bucket[bk]; j++)
                dT_tuple_count_after = dT_tuple_count_after + delta_sub_bucket_size[bk][j];
        }
        MPI_Allreduce(&dT_tuple_count_after, &dT_tuple_count_after_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        //if (rank == 0)
        //    std::cout << "[dT] After Load balancing " << dT_tuple_count_before_global << " " << dT_tuple_count_after_global << std::endl;
        assert(dT_tuple_count_before_global == dT_tuple_count_after_global);
#endif
        return true;
    }


    bool load_balance_merge_full(relation* rel, float rf)
    {
        u32 rank = (u32)mcomm.get_rank();
        int nprocs = mcomm.get_nprocs();
        int buckets = mcomm.get_number_of_buckets();
        u32* bucket_map = rel->get_bucket_map();
        u32* sub_bucket_count = rel->get_sub_bucket_count();
        google_relation* full = rel->get_full();
        u32** full_sub_bucket_size = rel->get_full_sub_bucket_element_count();
        u32** sub_bucket_rank = rel->get_sub_bucket_rank();
        int** distinct_sub_bucket_rank = rel->get_distinct_sub_bucket_rank();
        int* distinct_sub_bucket_rank_count = rel->get_distinct_sub_bucket_rank_count();
        int arity = rel->get_arity();

        u32 global_new_sub_bucket[buckets];
        memcpy(global_new_sub_bucket, sub_bucket_count, buckets * sizeof(u32));

        int mcount = 0;
        u32 new_sub_bucket_count = 0;
        u32 total_old_sub_bucket_count = 0;
        for (int i = 0; i < buckets; i++)
        {
            total_old_sub_bucket_count = total_old_sub_bucket_count + sub_bucket_count[i];
            if (sub_bucket_count[i] >= rf)
                mcount++;
        }


        if (mcount > 0.6* buckets)
        {
            for (int i = 0; i < buckets; i++)
            {
                global_new_sub_bucket[i] = global_new_sub_bucket[i] / rf;
                if (global_new_sub_bucket[i] == 0)
                    global_new_sub_bucket[i] = 1;

                new_sub_bucket_count = new_sub_bucket_count + global_new_sub_bucket[i];
            }
            if (rank == 0)
                std::cout << "[YES] Bucket Consolidation [" << mcount << " " << 0.6 * buckets << "] Old Sub Bucket count " << total_old_sub_bucket_count << " New Sub Bucket count " << new_sub_bucket_count << std::endl;
        }
        else
        {
            if (rank == 0)
                std::cout << "[NO] Bucket Consolidation [" << mcount << " " << 0.6 * buckets << "] Old Sub Bucket count " << total_old_sub_bucket_count << std::endl;
            return false;
        }

        int rcount = sub_bucket_rank[buckets-1][sub_bucket_count[buckets-1] - 1] + 1;
        for (int b = 0; b < buckets; b++)
        {
            if (sub_bucket_count[b] > global_new_sub_bucket[b])
            {
                bucket_map[b] = 0;

                std::unordered_set<int> distinct_t_ranks;
                int* temp = new int[sub_bucket_count[b]];
                memcpy(temp, sub_bucket_rank[b], sizeof(int) * sub_bucket_count[b]);

                delete[] sub_bucket_rank[b];
                sub_bucket_rank[b] = new u32[global_new_sub_bucket[b]];

                memcpy(sub_bucket_rank[b], temp, sizeof(int) * global_new_sub_bucket[b]);
                delete[] temp;

                for (u64 x = 0; x < global_new_sub_bucket[b]; x++)
                {
                    if (sub_bucket_rank[b][x] == rank)
                        bucket_map[b] = 1;

                    distinct_t_ranks.insert(sub_bucket_rank[b][x]);
                    rcount++;
                }

                delete[] distinct_sub_bucket_rank[b];
                distinct_sub_bucket_rank[b] = new int[distinct_t_ranks.size()];
                u32 x = 0;
                for ( auto it = distinct_t_ranks.begin(); it != distinct_t_ranks.end(); ++it )
                {
                    distinct_sub_bucket_rank[b][x] = *it;
                    x++;
                }
                distinct_sub_bucket_rank_count[b] = x;
            }

            //if (rank == 0)
            //    std::cout << "GNSB " << b << " " << global_new_sub_bucket[b] << std::endl;
        }


        u32* old_sub_bucket_count = new u32[buckets];
        memcpy(old_sub_bucket_count, sub_bucket_count, buckets * sizeof(u32));
        memcpy(sub_bucket_count, global_new_sub_bucket, sizeof(u32) * buckets);


        int process_size[nprocs];
        memset(process_size, 0, nprocs * sizeof(int));

        vector_buffer* process_data_vector = (vector_buffer*)malloc(sizeof(vector_buffer) * nprocs);
        for (int i = 0; i < nprocs; ++i) {
            process_data_vector[i] = vector_buffer_create_empty();
        }

        int prefix_sum_process_size[nprocs];
        memset(prefix_sum_process_size, 0, nprocs * sizeof(int));

        int recv_process_size_buffer[nprocs];
        memset(recv_process_size_buffer, 0, nprocs * sizeof(int));

        int process_size_dt[nprocs];
        memset(process_size_dt, 0, nprocs * sizeof(int));

        vector_buffer* process_data_vector_dt1 = (vector_buffer*)malloc(sizeof(vector_buffer) * nprocs);
        for (int i = 0; i < nprocs; ++i) {
            process_data_vector_dt1[i] = vector_buffer_create_empty();
        }

        int prefix_sum_process_size_dt[nprocs];
        memset(prefix_sum_process_size_dt, 0, nprocs * sizeof(int));

        int recv_process_size_buffer_dt[nprocs];
        memset(recv_process_size_buffer_dt, 0, nprocs * sizeof(int));

        u64 val[2] = {0, 0};
        for (int bk = 0; bk < buckets; bk++)
        {
            if (old_sub_bucket_count[bk] != global_new_sub_bucket[bk])
            {
                delete[] full_sub_bucket_size[bk];
                full_sub_bucket_size[bk] = new u32[global_new_sub_bucket[bk]];
                memset(full_sub_bucket_size[bk], 0, sizeof(u32) * global_new_sub_bucket[bk]);


                for (auto it = full[bk].begin(); it != full[bk].end(); it++)
                {
                    Map0* it2 = it->second;
                    for (auto dit2 = it2->begin(); dit2 != it2->end(); dit2++)
                    {
                        val[0] = it->first;
                        val[1] = dit2->first;
                        uint64_t bucket_id = hash_function(val[0]) % buckets;
                        uint64_t sub_bucket_id = hash_function(val[1]) % global_new_sub_bucket[bucket_id];

                        int index = sub_bucket_rank[bucket_id][sub_bucket_id];
                        process_size[index] = process_size[index] + arity;


                        vector_buffer_append(&process_data_vector[index],(unsigned char *) val, sizeof(u64) * arity);

                        //if (index == 0)
                        //std::cout << "BI and SBI " << val[0] << " " << val[1] << " " << hash_function(val[0]) << " " << bucket_id << " " << sub_bucket_id << " rank " << index << std::endl;
                    }
                }

                for(google_relation::iterator iy2 = full[bk].begin(); iy2 != full[bk].end(); iy2++)
                    delete (iy2->second);

                full[bk].clear();
            }
        }

        u64 outer_hash_buffer_size = 0;
        u64* outer_hash_data;
        all_to_all_comm(process_data_vector, process_size, &outer_hash_buffer_size, &outer_hash_data);

        free (process_data_vector);

        //if (rank == 0)
        //    std::cout << "outer_hash_buffer_size " << outer_hash_buffer_size <<std::endl;

        u64 t[2];
        for (u64 in = 0; in < outer_hash_buffer_size; in = in + arity)
        {
            t[0] = outer_hash_data[in];
            t[1] = outer_hash_data[in + 1];
            rel->insert_in_full(t);
        }
        delete[] outer_hash_data;
        delete[] old_sub_bucket_count;

        return true;
    }


    void set_refinement_factor(double rf)
    {
        refinement_factor = rf;
    }

    void set_refinement_interval(int ri)
    {
        refinement_ts = ri;
    }


    void set_threshold(int thold)
    {
        threshold = thold;
    }

    void enable_comm_compaction()
    {
        comm_compaction = true;
    }

    void disable_comm_compaction()
    {
        comm_compaction = false;
    }

    void execute()
    {
        //u32 iteration = 0;

        //u32 clique_comm_count = 0;
        u32 local_join_count = 0;
        //u32 all_to_all_count = 0;
        //u32 insert_in_full_count = 0;
        //u32 insert_in_newt_count = 0;

        //u32 running_clique_comm_count = 0;
        //u32 running_local_join_count = 0;
        //u32 running_all_to_all_count = 0;
        //u32 running_insert_in_full_count = 0;
        //u32 running_insert_in_newt_count = 0;

        double running_clique_comm_time = 0;
        double running_local_join_time = 0;
        double running_all_to_all_time = 0;
        double running_insert_in_full_time = 0;
        double running_insert_in_newt_time = 0;
        double running_verify_time = 0;
        double running_lb = 0;

        double clique_start = 0, clique_end = 0;
        double local_join_start = 0, local_join_end = 0;
        double all_to_all_start = 0, all_to_all_end = 0;
        double insert_full_start = 0, insert_full_end = 0;
        double insert_newt_start = 0, insert_newt_end = 0;
        double verify_start = 0, verify_end = 0;
        double lb_start = 0, lb_end = 0;

        double iteration_time = 0;
        double running_time = 0;

        u32 RA_count = RA_list.size();
        int *offset = new int[RA_count];
        for (u32 i =0; i < RA_count; i++)
            offset[i] = 0;

        bool local_join_status = true;
        int rank = mcomm.get_rank();
        int iteration = 1;

        double start_time = MPI_Wtime();

        int outer_loop = 0;
        int inner_loop = 0;



        int itx = 10;
        if (comm_compaction== false)
        {
            if (rank == 0)
                std::cout << "Initial Load balancing " << std::endl;
            initial_load_balance(refinement_factor);
            if (rank == 0)
                std::cout << std::endl << std::endl;

            //while (itx != 0)
            while (true)
            {
                clique_start = MPI_Wtime();
                clique_comm(local_join_status);
                clique_end = MPI_Wtime();
                running_clique_comm_time = running_clique_comm_time + (clique_end - clique_start);

                local_join_start = MPI_Wtime();
                local_join_status = local_join(&local_join_count, offset, iteration);
                local_join_end = MPI_Wtime();
                running_local_join_time = running_local_join_time + (local_join_end - local_join_start);

                all_to_all_start = MPI_Wtime();
                all_to_all();
                all_to_all_end = MPI_Wtime();
                running_all_to_all_time = running_all_to_all_time + (all_to_all_end - all_to_all_start);

                insert_newt_start = MPI_Wtime();
                local_insert_in_newt();
                insert_newt_end = MPI_Wtime();
                running_insert_in_newt_time = running_insert_in_newt_time + (insert_newt_end - insert_newt_start);


                insert_full_start = MPI_Wtime();
                local_insert_in_full(local_join_status);
                insert_full_end = MPI_Wtime();
                running_insert_in_full_time = running_insert_in_full_time + (insert_full_end - insert_full_start);


                lb_start = MPI_Wtime();
                if (refinement_ts != 0)
                    if (outer_loop % refinement_ts == 0)
                        load_balance(local_join_status, refinement_factor);

                if (local_join_status == true)
                {
                    inner_loop = 0;
                    outer_loop++;
                }
                lb_end = MPI_Wtime();
                running_lb = running_lb + (lb_end - lb_start);


                verify_start = MPI_Wtime();
                if (outer_loop % 10 == 0)
                {
                    print_full();
                    if (check_for_fixed_point(local_join_status) == true)
                    {
                        verify_end = MPI_Wtime();
                        running_verify_time = running_verify_time + (verify_end - verify_start);
                        iteration_time = (verify_end - verify_start) + (lb_end - lb_start) + (insert_newt_end - insert_newt_start) + (insert_full_end - insert_full_start) + (all_to_all_end - all_to_all_start) + (local_join_end - local_join_start) + (clique_end - clique_start);
                        running_time = running_time + iteration_time;

                        if (rank == 0)
                            std::cout << "T NCC OL " << outer_loop << " IL " << inner_loop << " [" << iteration << "] TH " << threshold << " RF " << refinement_factor << " RI " << refinement_ts << " " << running_time << " " << iteration_time
                                      << " Clique " <<  (clique_end - clique_start)
                                      << " Local Join " <<  (local_join_end - local_join_start)
                                      << " A2A " <<  (all_to_all_end - all_to_all_start)
                                      << " Insert Full " <<  (insert_full_end - insert_full_start)
                                      << " Insert newt " <<  (insert_newt_end - insert_newt_start)
                                      << " LB " <<  (lb_end - lb_start)
                                      << " Verify " <<  (verify_end - verify_start)
                                      << std::endl
                                      << std::endl;
                        break;
                    }
                }
                verify_end = MPI_Wtime();
                running_verify_time = running_verify_time + (verify_end - verify_start);

                iteration_time = (verify_end - verify_start) + (lb_end - lb_start) + (insert_newt_end - insert_newt_start) + (insert_full_end - insert_full_start) + (all_to_all_end - all_to_all_start) + (local_join_end - local_join_start) + (clique_end - clique_start);
                running_time = running_time + iteration_time;

                if (rank == 0)
                    std::cout << "F NCC OL " << outer_loop << " IL " << inner_loop << " [" << iteration << "] TH " << threshold << " RF " << refinement_factor << " RI " << refinement_ts << " " << running_time << " " << iteration_time
                              << " Clique " <<  (clique_end - clique_start)
                              << " Local Join " <<  (local_join_end - local_join_start)
                              << " A2A " <<  (all_to_all_end - all_to_all_start)
                              << " Insert Full " <<  (insert_full_end - insert_full_start)
                              << " Insert newt " <<  (insert_newt_end - insert_newt_start)
                              << " LB " <<  (lb_end - lb_start)
                              << " Verify " <<  (verify_end - verify_start)
                              << std::endl
                              << std::endl
                              << std::endl;

                inner_loop++;
                iteration++;
                itx--;
            }
        }
#if 0
        else
        {
            while (true)
            {
                clique_start = MPI_Wtime();
                clique_comm(local_join_status);
                clique_end = MPI_Wtime();
                running_clique_comm_time = running_clique_comm_time + (clique_end - clique_start);

                local_join_start = MPI_Wtime();
                local_join_status = cumulative_local_join(&local_join_count, offset, iteration);
                local_join_end = MPI_Wtime();
                running_local_join_time = running_local_join_time + (local_join_end - local_join_start);


                all_to_all_start = MPI_Wtime();
                cumulative_all_to_all();
                all_to_all_end = MPI_Wtime();
                running_all_to_all_time = running_all_to_all_time + (all_to_all_end - all_to_all_start);


                insert_newt_start = MPI_Wtime();
                cumulative_local_insert_in_newt();
                insert_newt_end = MPI_Wtime();
                running_insert_in_newt_time = running_insert_in_newt_time + (insert_newt_end - insert_newt_start);


                insert_full_start = MPI_Wtime();
                local_insert_in_full(local_join_status);
                insert_full_end = MPI_Wtime();
                running_insert_in_full_time = running_insert_in_full_time + (insert_full_end - insert_full_start);

                if (local_join_status == true)
                {
                    inner_loop = 0;
                    outer_loop++;
                }

                //if (rank == 0)
                //std::cout << "ITERATION [" << iteration <<"] " <<  local_join_count << " " << all_to_all_count << " " << insert_in_full_count << std::endl;
                //std::cout << "RUNNING ITERATION [" << iteration <<"] " << running_clique_comm_count << " " << running_local_join_count << " " << running_all_to_all_count << " " << running_insert_in_full_count << " " << running_insert_in_delta << std::endl;


                verify_start = MPI_Wtime();
                if (check_for_fixed_point(local_join_status) == true)
                {
                    verify_end = MPI_Wtime();
                    running_verify_time = running_verify_time + (verify_end - verify_start);
                    iteration_time = (verify_end - verify_start) + (insert_newt_end - insert_newt_start) + (insert_full_end - insert_full_start) + (all_to_all_end - all_to_all_start) + (local_join_end - local_join_start) + (clique_end - clique_start);
                    running_time = running_time + iteration_time;

                    if (rank == 0)
                        std::cout << "T CC OL " << outer_loop << " IL " << inner_loop << " [" << iteration << "] " << running_time << " " << iteration_time
                                  << " Clique " <<  (clique_end - clique_start)
                                  << " LJ " <<  (local_join_end - local_join_start)
                                  << " All to All " <<  (all_to_all_end - all_to_all_start)
                                  << " Insert Full " <<  (insert_full_end - insert_full_start)
                                  << " Insert newt " <<  (insert_newt_end - insert_newt_start)
                                  << " Verify " <<  (verify_end - verify_start)
                                  << std::endl
                                  << std::endl;
                    break;
                }
                verify_end = MPI_Wtime();
                running_verify_time = running_verify_time + (verify_end - verify_start);

                iteration_time = (verify_end - verify_start) + (insert_newt_end - insert_newt_start) + (insert_full_end - insert_full_start) + (all_to_all_end - all_to_all_start) + (local_join_end - local_join_start) + (clique_end - clique_start);
                running_time = running_time + iteration_time;

                if (rank == 0)
                    std::cout << "F CC OL " << outer_loop << " IL " << inner_loop << " [" << iteration << "] " << running_time << " " << iteration_time
                              << " Clique " <<  (clique_end - clique_start)
                              << " LJ " <<  (local_join_end - local_join_start)
                              << " All to All " <<  (all_to_all_end - all_to_all_start)
                              << " Insert Full " <<  (insert_full_end - insert_full_start)
                              << " Insert newt " <<  (insert_newt_end - insert_newt_start)
                              << " Verify " <<  (verify_end - verify_start)
                              << std::endl
                              << std::endl
                              << std::endl;

                inner_loop++;
                iteration++;
            }
        }
#endif
        double end_time = MPI_Wtime();
        delete[] offset;

        if (rank == 0)
            std::cout << "Threshold " << threshold
                      << " Total Time: [" << (end_time - start_time)
                      << " " << running_time
                      << " "
                      << (running_clique_comm_time + running_local_join_time + running_all_to_all_time + running_insert_in_newt_time + running_insert_in_full_time + running_lb + running_verify_time)
                      << "] Clique " << running_clique_comm_time
                      << " LJ " << running_local_join_time
                      << " A2A " << running_all_to_all_time
                      << " Insert in new " << running_insert_in_newt_time
                      << " Insert in full " << running_insert_in_full_time
                      << " LB " << running_lb
                      << " FPC " << running_verify_time << " ";
        print_full();
    }

};

#endif
