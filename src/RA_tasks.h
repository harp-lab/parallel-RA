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


class RAM
{
private:
    bool comm_compaction;
    u32 threshold;

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

            if (mcomm.get_rank() == 0)
                std::cout << "ALWAYS " << std::endl;

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
                break;
            }
        }

        return local_insert;
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

        double clique_start = 0, clique_end = 0;
        double local_join_start = 0, local_join_end = 0;
        double all_to_all_start = 0, all_to_all_end = 0;
        double insert_full_start = 0, insert_full_end = 0;
        double insert_newt_start = 0, insert_newt_end = 0;
        double verify_start = 0, verify_end = 0;

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



        if (comm_compaction== false)
        {
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
                        std::cout << "T NCC OL " << outer_loop << " IL " << inner_loop << " [" << iteration << "] " << running_time << " " << iteration_time
                                  << " Clique " <<  (clique_end - clique_start)
                                  << " Local Join " <<  (local_join_end - local_join_start)
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
                    std::cout << "F NCC OL " << outer_loop << " IL " << inner_loop << " [" << iteration << "] " << running_time << " " << iteration_time
                              << " Clique " <<  (clique_end - clique_start)
                              << " Local Join " <<  (local_join_end - local_join_start)
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
                                  << " Local Join " <<  (local_join_end - local_join_start)
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
                              << " Local Join " <<  (local_join_end - local_join_start)
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
        double end_time = MPI_Wtime();
        delete[] offset;

        if (rank == 0)
            std::cout << "Threshold " << threshold << " Total Time: [" << (end_time - start_time) << " " << running_time << " " << (running_clique_comm_time + running_local_join_time + running_all_to_all_time + running_insert_in_newt_time + running_insert_in_full_time + running_verify_time) << "] Clique " << running_clique_comm_time << " Local Join " << running_local_join_time << " All to All " << running_all_to_all_time << " Insert in new " << running_insert_in_newt_time << " Insert in full " << running_insert_in_full_time << " Fixed point check " << running_verify_time << " ";
        print_full();
    }

};

#endif
