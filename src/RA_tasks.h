#ifndef RAM_TASKS_H
#define RAM_TASKS_H

#include <mpi.h>
#include <vector>
#include <unordered_set>
#include "balanced_hash_relation.h"
#include "balanced_parallel_io.h"
#include "parallel_join.h"
#include "btree/btree_map.h"
#include <queue>

u32 running_total = 0;
bool threshold_enabled = false;


class RAM
{
private:
    bool comm_compaction;
    u32 threshold;

    int refinement_chooser;
    double refinement_factor;

    u32 refinement_ts;

    std::vector<relation*> relation_manager;
    std::vector<parallel_RA*> RA_list;

    u64 *intra_bucket_buf_output_size;
    u64 **intra_bucket_buf_output;

    // local_compute_output is of size RA_list size x nprocs
    // this contains data per RA rule for per process
    vector_buffer** local_compute_output;

    // local_compute_output_size is of size RA_list size x nprocs
    // this contains size of data per RA rule for per process
    int **local_compute_output_size;

    // cumulative_tuple_process_map is of size nprocs
    // cumulative_tuple_process_map[i] contains number of tuples to be transmitted to rank-i process, across all Relation Algebra (RA) rules
    int *cumulative_tuple_process_map;

    u64 *cumulative_all_to_all_buffer;
    int* cumulative_all_to_all_recv_process_size_array;

    mpi_comm mcomm;


public:

    RAM()
    {
        this->refinement_ts = 10;
        this->comm_compaction = false;
        this->threshold = 1000000;
        this->refinement_factor = 4;
        this->refinement_chooser = 0;
    }

    void set_comm(mpi_comm mcomm)   {this->mcomm = mcomm;}
    void push_back(parallel_RA* pj) {RA_list.push_back(pj);}
    void push_relation(relation* G) {relation_manager.push_back(G);}
    void set_refinement_chooser(int rc) {refinement_chooser = rc;}
    void set_refinement_factor(double rf)   {refinement_factor = rf;}
    void set_refinement_interval(int ri)    {refinement_ts = ri;}
    void set_threshold(int thold)   {threshold = thold;}
    void enable_comm_compaction()   {comm_compaction = true;}
    void disable_comm_compaction()  {comm_compaction = false;}



    // Print relation size
    // this function is hacky right now
    void print_full()
    {
        u32 counter = 0;
        int rank = mcomm.get_rank();

        for (std::vector<relation*>::iterator it = relation_manager.begin() ; it != relation_manager.end(); ++it)
        {
            relation* rel = *it;
            u64 send = rel->get_full_element_count();
            u64 recv;
            MPI_Allreduce(&send, &recv, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
            if (mcomm.get_rank() == 0)
                std::cout << "[" << rank << "] Relation " << counter << " Full [" << recv << "] " << send << " Delta " << rel->get_delta_element_count() << std::endl;
            counter++;

#if 0
            if (rel->get_filename() != NULL)
            {
                std::vector<u64> prefix = {};
                vector_buffer vb = vector_buffer_create_empty();
                google_relation* full = rel->get_full();

                full[mcomm.get_rank()].as_vector_buffer(&vb, prefix, 0);

                if (mcomm.get_rank() == 0)
                {
                    mkdir(rel->get_filename(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

                    FILE* fp;
                    char metafname[PATH_MAX];
                    sprintf(metafname, "%s%s", rel->get_filename(), "/meta_data.txt");
                    //std::cout << "Metadta filename1 " << metafname << std::endl;
                    fp = fopen (metafname, "w");
                    fprintf (fp, "(row count)\n%d\n(col count)\n2",(int)recv);
                    fclose (fp);
                }

                MPI_Barrier(MPI_COMM_WORLD);

                int int_send = send;
                int scan;
                MPI_Scan(&int_send, &scan, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                MPI_File fp;
                MPI_Status status;
                MPI_Offset offset = scan;
                char fname[PATH_MAX];
                sprintf(fname, "%s%s", rel->get_filename(), "/data.raw");
                std::cout << "fname " << fname << "Offset " << offset << " send " << int_send << std::endl;
                MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
                MPI_File_write_at_all(fp, offset, vb.buffer, int_send, MPI_UNSIGNED_LONG_LONG, &status);
                MPI_File_close(&fp);
            }
#endif

        }
        return;
    }


    // check for fixed point
    bool check_for_fixed_point()
    {
        bool fixed_point = true;
        for (std::vector<relation*>::iterator it = relation_manager.begin() ; it != relation_manager.end(); ++it)
        {
            relation* rel = *it;
            fixed_point = fixed_point & rel->fixed_point_check();
        }

        return fixed_point;
    }


    // intra_bucket Comm
    u64 intra_bucket_comm()
    {
        u64 total_data_moved = 0;
        u32 counter = 0;
        u32 RA_count = RA_list.size();

        intra_bucket_buf_output_size = new u64[RA_count];
        intra_bucket_buf_output = new u64*[RA_count];


        for (std::vector<parallel_RA*>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
        {
            parallel_RA* current_ra = *it;

            // No intra-bucket comm required for copy
            if (current_ra->get_RA_type() == COPY)
                continue;


            // No intra-bucket comm required for acopy
            if (current_ra->get_RA_type() == ACOPY)
                continue;

            // Intra-bucket comm for joins
            else if (current_ra->get_RA_type() == JOIN)
            {
                // Join between delta and delta
                if (current_ra->get_join_input0_graph_type() == DELTA && current_ra->get_join_input1_graph_type() == DELTA)
                {
                    // We are picking the second delta as the source (ROOM FOR OPTIMIZATION HERE)
                    relation* input = current_ra->get_join_input1();
                    relation* output = current_ra->get_join_input0();

                    current_ra->intra_bucket_comm(input->get_delta(),
                                                  input->get_number_of_buckets(), input->get_distinct_sub_bucket_rank_count(), input->get_distinct_sub_bucket_rank(), input->get_bucket_map(),
                                                  output->get_number_of_buckets(), output->get_distinct_sub_bucket_rank_count(), output->get_distinct_sub_bucket_rank(), output->get_bucket_map(),
                                                  &intra_bucket_buf_output_size[counter], &intra_bucket_buf_output[counter]);
                    total_data_moved = total_data_moved + intra_bucket_buf_output_size[counter];
                }

                // Join between delta and full
                else if (current_ra->get_join_input0_graph_type() == DELTA && current_ra->get_join_input1_graph_type() == FULL)
                {
                    relation* input = current_ra->get_join_input0();
                    relation* output = current_ra->get_join_input1();

                    current_ra->intra_bucket_comm(input->get_delta(),
                                                  input->get_number_of_buckets(), input->get_distinct_sub_bucket_rank_count(), input->get_distinct_sub_bucket_rank(), input->get_bucket_map(),
                                                  output->get_number_of_buckets(), output->get_distinct_sub_bucket_rank_count(), output->get_distinct_sub_bucket_rank(), output->get_bucket_map(),
                                                  &intra_bucket_buf_output_size[counter], &intra_bucket_buf_output[counter]);
                    total_data_moved = total_data_moved + intra_bucket_buf_output_size[counter];

                }

                // Join between full and delta
                else if (current_ra->get_join_input0_graph_type() == FULL && current_ra->get_join_input1_graph_type() == DELTA)
                {
                    relation* input = current_ra->get_join_input1();
                    relation* output = current_ra->get_join_input0();

                    current_ra->intra_bucket_comm(input->get_delta(),
                                                  input->get_number_of_buckets(), input->get_distinct_sub_bucket_rank_count(), input->get_distinct_sub_bucket_rank(), input->get_bucket_map(),
                                                  output->get_number_of_buckets(), output->get_distinct_sub_bucket_rank_count(), output->get_distinct_sub_bucket_rank(), output->get_bucket_map(),
                                                  &intra_bucket_buf_output_size[counter], &intra_bucket_buf_output[counter]);
                    total_data_moved = total_data_moved + intra_bucket_buf_output_size[counter];
                }

                // Join between full and full
                else if (current_ra->get_join_input0_graph_type() == FULL && current_ra->get_join_input1_graph_type() == FULL)
                {
                    relation* input = current_ra->get_join_input1();
                    relation* output = current_ra->get_join_input0();

                    current_ra->intra_bucket_comm(input->get_full(),
                                                  input->get_number_of_buckets(), input->get_distinct_sub_bucket_rank_count(), input->get_distinct_sub_bucket_rank(), input->get_bucket_map(),
                                                  output->get_number_of_buckets(), output->get_distinct_sub_bucket_rank_count(), output->get_distinct_sub_bucket_rank(), output->get_bucket_map(),
                                                  &intra_bucket_buf_output_size[counter], &intra_bucket_buf_output[counter]);
                    total_data_moved = total_data_moved + intra_bucket_buf_output_size[counter];
                }
            }
            counter++;
        }

        return total_data_moved;
    }


    bool local_compute(u32* local_join_count, int* offset)
    {
        u32 join_tuples = 0;
        u32 total_join_tuples = 0;
        u32 RA_count = RA_list.size();
        u32 nprocs = (u32)mcomm.get_nprocs();

        bool join_completed = true;

        local_compute_output = new vector_buffer*[RA_count];
        local_compute_output_size = new int*[RA_count];

        cumulative_tuple_process_map = new int[nprocs];
        memset(cumulative_tuple_process_map, 0, nprocs * sizeof(int));

        for (u32 i = 0; i < RA_count; i++)
        {
            local_compute_output[i] = new vector_buffer[nprocs];
            for (u32 j = 0; j < nprocs; ++j) {
                local_compute_output[i][j] = vector_buffer_create_empty();
            }

            local_compute_output_size[i] = new int[nprocs];
            memset(local_compute_output_size[i], 0, nprocs * sizeof(int));
        }

        u32 counter = 0;
        for (std::vector<parallel_RA*>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
        {
            if ((*it)->get_RA_type() == COPY)
            {
                parallel_copy* current_ra = (parallel_copy*) *it;


                int reorder_map_array_size;
                int* reorder_map_array;
                current_ra->get_copy_rename_index(&reorder_map_array, &reorder_map_array_size);

                if (current_ra->get_copy_input0_graph_type() == DELTA)
                {
                    relation* output_relation = current_ra->get_copy_output();
                    relation* input_relation = current_ra->get_copy_input();

                    current_ra->local_copy(input_relation->get_delta(), input_relation->get_bucket_map(),
                                           output_relation,
                                           reorder_map_array_size, reorder_map_array,
                                           local_compute_output[counter],
                                           local_compute_output_size[counter],
                                           cumulative_tuple_process_map);
                }
                if (current_ra->get_copy_input0_graph_type() == FULL)
                {
                    relation* output_relation = current_ra->get_copy_output();
                    relation* input_relation = current_ra->get_copy_input();

                    current_ra->local_copy(input_relation->get_full(), input_relation->get_bucket_map(),
                                           output_relation,
                                           reorder_map_array_size, reorder_map_array,
                                           local_compute_output[counter],
                                           local_compute_output_size[counter],
                                           cumulative_tuple_process_map);
                }
            }

            if ((*it)->get_RA_type() == ACOPY)
            {

            }

            else if ((*it)->get_RA_type() == JOIN)
            {
                parallel_join* current_ra = (parallel_join*) *it;
                relation* output_relation = current_ra->get_join_output();
                int join_column_count = current_ra->get_join_column_count();

                int reorder_map_array_size;
                int* reorder_map_array;
                current_ra->get_join_projection_index(&reorder_map_array, &reorder_map_array_size);

                if (current_ra->get_join_input0_graph_type() == DELTA && current_ra->get_join_input1_graph_type() == DELTA)
                {
                    relation* delta0 = current_ra->get_join_input0();
                    relation* delta1 = current_ra->get_join_input1();

                    join_completed = join_completed &
                            current_ra->local_join(intra_bucket_buf_output_size[counter], delta1->get_arity(), intra_bucket_buf_output[counter], 0,
                                                   delta0->get_delta(), delta0->get_delta_element_count(), delta0->get_arity(),
                                                   reorder_map_array_size, reorder_map_array,
                                                   output_relation,
                                                   local_compute_output, local_compute_output_size, cumulative_tuple_process_map,
                                                   threshold, &(offset[counter]), join_column_count,
                                                   &join_tuples, counter);
                    total_join_tuples = total_join_tuples + join_tuples;
                }
                else if (current_ra->get_join_input0_graph_type() == DELTA && current_ra->get_join_input1_graph_type() == FULL)
                {
                    relation* full = current_ra->get_join_input1();
                    relation* delta = current_ra->get_join_input0();


                    join_completed = join_completed &
                            current_ra->local_join(intra_bucket_buf_output_size[counter], delta->get_arity(), intra_bucket_buf_output[counter], 0,
                                                   full->get_full(), full->get_full_element_count(), full->get_arity(),
                                                   reorder_map_array_size, reorder_map_array,
                                                   output_relation,
                                                   local_compute_output, local_compute_output_size, cumulative_tuple_process_map,
                                                   threshold,  &(offset[counter]), join_column_count,
                                                   &join_tuples, counter);
                    total_join_tuples = total_join_tuples + join_tuples;
                }
                else if (current_ra->get_join_input0_graph_type() == FULL && current_ra->get_join_input1_graph_type() == DELTA)
                {
                    relation* full = current_ra->get_join_input0();
                    relation* delta = current_ra->get_join_input1();

                    join_completed = join_completed &
                            current_ra->local_join(intra_bucket_buf_output_size[counter], delta->get_arity(), intra_bucket_buf_output[counter], 1,
                                                   full->get_full(), full->get_full_element_count(), full->get_arity(),
                                                   reorder_map_array_size, reorder_map_array,
                                                   output_relation,
                                                   local_compute_output, local_compute_output_size, cumulative_tuple_process_map,
                                                   threshold,  &(offset[counter]), join_column_count,
                                                   &join_tuples, counter);
                    total_join_tuples = total_join_tuples + join_tuples;
                }
            }
            counter++;
        }

        /*

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
                    delete[] intra_bucket_buf_output[counter];

                offset[counter] = 0;
                counter++;
            }

            delete[] intra_bucket_buf_output_size;
            delete[] intra_bucket_buf_output;
            return true;
        }
        else
            return false;
        */

        counter = 0;
        for (std::vector<parallel_RA*>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
        {
            parallel_RA* current_ra = *it;
            if (current_ra->get_RA_type() == JOIN)
                delete[] intra_bucket_buf_output[counter];
        }

        delete[] intra_bucket_buf_output_size;
        delete[] intra_bucket_buf_output;

        return true;
    }



    int all_to_all()
    {
        int counter = 0;
        u64 total_all_to_all = 0;
        u32 RA_count = RA_list.size();
        int nprocs = mcomm.get_nprocs();

        cumulative_all_to_all_recv_process_size_array = new int[RA_count * nprocs];

        int *flat_process_size = new int[nprocs * RA_count];
        for (int i = 0; i < nprocs; i++)
            for (u32 j = 0; j < RA_count; j++)
                flat_process_size[counter++] = local_compute_output_size[j][i];

        /* This step prepares for actual data transfer */
        /* Every process sends to every other process the amount of data it is going to send */
        //int recv_process_size_array[nprocs * RA_count];
        memset(cumulative_all_to_all_recv_process_size_array, 0, RA_count * nprocs * sizeof(int));
        MPI_Alltoall(flat_process_size, RA_count, MPI_INT, cumulative_all_to_all_recv_process_size_array, RA_count, MPI_INT, mcomm.get_comm());

        delete[] flat_process_size;


        // PREPARING THE SEND PREFIX
        int cumulative_prefix_sum_process_size[nprocs];
        memset(cumulative_prefix_sum_process_size, 0, nprocs * sizeof(int));

        for (int i = 1; i < nprocs; i++)
            cumulative_prefix_sum_process_size[i] = cumulative_prefix_sum_process_size[i - 1] + cumulative_tuple_process_map[i - 1];

        int process_data_buffer_size = cumulative_prefix_sum_process_size[nprocs - 1] + cumulative_tuple_process_map[nprocs - 1];


        // PREPARING THE SEND DATA
        u64* process_data = 0;
        process_data = new u64[process_data_buffer_size];
        memset(process_data, 0, process_data_buffer_size * sizeof(u64));

        u32 boffset = 0;
        for(int i = 0; i < nprocs; i++)
        {
            for (u32 r = 0; r < RA_count; r++)
            {
                memcpy(process_data + boffset, (&local_compute_output[r][i])->buffer, (&local_compute_output[r][i])->size);
                boffset = boffset + ((&local_compute_output[r][i])->size)/8;
                vector_buffer_free(&local_compute_output[r][i]);
            }
        }

        // PREPARING THE RECEIVE PREFIX
        int cumulative_recv_process_size_array[nprocs];
        memset(cumulative_recv_process_size_array, 0, nprocs * sizeof(int));

        for (u32 k = 0; k < RA_count * nprocs; k++)
            cumulative_recv_process_size_array[k/RA_count] = cumulative_recv_process_size_array[k/RA_count] + cumulative_all_to_all_recv_process_size_array[k];

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
        cumulative_all_to_all_buffer = new u64[outer_hash_buffer_size];
        memset(cumulative_all_to_all_buffer, 0, outer_hash_buffer_size * sizeof(u64));


        // ALL TO ALL COMM
        MPI_Alltoallv(process_data, cumulative_tuple_process_map, cumulative_prefix_sum_process_size, MPI_UNSIGNED_LONG_LONG, cumulative_all_to_all_buffer, cumulative_recv_process_size_array, cumulative_prefix_sum_recv_process_size_buffer, MPI_UNSIGNED_LONG_LONG, mcomm.get_comm());

        delete[] process_data;

        /*
        if (mcomm.get_rank() == 1)
        {
            for (int x; x < outer_hash_buffer_size; x=x+2)
                std::cout << "XXXXXXXXXX " << cumulative_all_to_all_buffer[x] << " " << cumulative_all_to_all_buffer[x+1] << std::endl;
        }
        */

        for (u32 i = 0; i < RA_count; i++)
        {
            delete[] local_compute_output[i];
            delete[] local_compute_output_size[i];
        }

        delete[] local_compute_output;
        delete[] local_compute_output_size;
        delete[] cumulative_tuple_process_map;

        return total_all_to_all;
    }


    void local_insert_in_newt()
    {
        u32 successful_insert = 0;
        int nprocs = mcomm.get_nprocs();
        //int offssizeet = 0;
        int RA_count = RA_list.size();


        // RA_count = 4
        // nprocs = 16
        // RA_count * nprocs = 64
        // 10 10 10 10
        // k = 0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60

        u32 starting = 0;
        for (int k = 0; k < RA_count * nprocs; k++)
        {
            u32 ra_id = k % RA_count;
            u32 elements_to_read = cumulative_all_to_all_recv_process_size_array[k];
            relation* output;

            if (RA_list[ra_id]->get_RA_type() == COPY)
                output = RA_list[ra_id]->get_copy_output();
            else
                output = RA_list[ra_id]->get_join_output();

            u32 arity = output->get_arity();

            //if (mcomm.get_rank() == 1)
            //    std::cout << "elements_to_read " << elements_to_read << std::endl;

            for (u32 x = starting; x < starting + elements_to_read; x=x+arity)
            {
                u64 t[arity];
                for (u32 a = 0; a < arity; a++)
                    t[a] = cumulative_all_to_all_buffer[x + a];

                //if (mcomm.get_rank() == 1)
                //    std::cout << "[Trying] to insert in newt " << t[0] << " " << t[1] << std::endl;

                if (output->find_in_full(t) == false)
                {
                    if (output->insert_in_newt(t) == true)
                    {
                        //if (mcomm.get_rank() == 1)
                        //    std::cout << "[Success] Insert in newt " << t[0] << " " << t[1] << std::endl;
                        successful_insert++;
                    }
                }
            }
            starting = starting + elements_to_read;

        }


        /*
        for (int k = 0; k < RA_count * nprocs; k= k + RA_count)
        {
            for (int i = 0; i < RA_count; i++)
            {
                if (RA_list[i]->get_RA_type() == COPY)
                {
                    relation* output = RA_list[i]->get_copy_output();
                    u32 arity = output->get_arity();

                    for (int j = 0; j < cumulative_all_to_all_recv_process_size_array[k + i]; j = j + arity)
                    {
                        u64 t[arity];
                        for (u32 a = 0; a < arity; a++)
                            t[a] = cumulative_all_to_all_buffer[j + a];

                        if (mcomm.get_rank() == 0)
                            std::cout << "[Trying] to insert in newt " << t[0] << " " << t[1] << std::endl;

                        if (output->find_in_full(t) == false)
                        {
                            if (output->insert_in_newt(t) == true)
                            {
                                if (mcomm.get_rank() == 0)
                                    std::cout << "[Success] Insert in newt " << t[0] << " " << t[1] << std::endl;
                                successful_insert++;
                            }
                        }
                    }
                }
                else if (RA_list[i]->get_RA_type() == JOIN)
                {
                    relation* output = RA_list[i]->get_join_output();
                    u32 arity = output->get_arity();

                    for (int j = 0; j < cumulative_all_to_all_recv_process_size_array[k + i]; j = j + arity)
                    {
                        u64 t[arity];
                        for (u32 a = 0; a < arity; a++)
                            t[a] = cumulative_all_to_all_buffer[j + a];

                        if (output->find_in_full(t) == false)
                        {
                            output->insert_in_newt(t);
                            successful_insert++;
                        }
                    }
                }

                offset = offset + cumulative_all_to_all_recv_process_size_array[k + i];
            }
        }
        */

        delete[] cumulative_all_to_all_recv_process_size_array;
        delete[] cumulative_all_to_all_buffer;
    }


    u32 local_insert_in_full()
    {
        u32 local_insert = 0;

        for (std::vector<relation*>::iterator it = relation_manager.begin() ; it != relation_manager.end(); ++it)
        {
            relation* current_r = *it;
            //current_r->print();
            current_r->insert_delta_in_full();
            current_r->local_insert_in_delta();
            //

        }

#if 0
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
                break;

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
#endif

        return local_insert;
    }


#if 0
    void load_balance(float r_factor, int rc)
    {
        //if (threshold_reached == false)
        //    return;

        // load balance all relations
        for (std::vector<relation*>::iterator it = relation_manager.begin() ; it != relation_manager.end(); ++it)
        {
            relation* current_relation = *it;
            if (load_balance_merge_full_and_delta(current_relation, r_factor) == false)
                load_balance_split_full_and_delta(current_relation, r_factor, rc);

            u64 max_full_element_count = 0;
            u64 min_full_element_count = 0;
            u64 sum_full_element_count = 0;
            u64 full_element_count = 0;

            u64 max_delta_element_count = 0;
            u64 min_delta_element_count = 0;
            u64 sum_delta_element_count = 0;
            u64 delta_element_count = 0;

            u32* bucket_map = current_relation->get_bucket_map();
            u32* sub_bucket_count = current_relation->get_sub_bucket_count();

            u32** full_sub_bucket_size = current_relation->get_full_sub_bucket_element_count();
            u32** delta_sub_bucket_size = current_relation->get_delta_sub_bucket_element_count();


            for (int i = 0; i < mcomm.get_number_of_buckets(); i++)
            {
                if (bucket_map[i] == 1)
                {
                    for (u32 j = 0; j < sub_bucket_count[i]; j++)
                    {
                        full_element_count = full_element_count + full_sub_bucket_size[i][j];
                        delta_element_count = delta_element_count + delta_sub_bucket_size[i][j];
                    }
                }
            }

            MPI_Allreduce(&full_element_count, &max_full_element_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(&full_element_count, &min_full_element_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(&full_element_count, &sum_full_element_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

            MPI_Allreduce(&delta_element_count, &max_delta_element_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(&delta_element_count, &min_delta_element_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(&delta_element_count, &sum_delta_element_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

            if (mcomm.get_rank() == 0)
            {
                if (min_full_element_count != 0 && min_delta_element_count != 0)
                    std::cout << "[LOAD BALANCING] Full (" << max_full_element_count << ", " << min_full_element_count << ", " << sum_full_element_count/mcomm.get_nprocs() << ") "
                              << (float) max_full_element_count/min_full_element_count << " "
                              << "Delta (" << max_delta_element_count << ", " << min_delta_element_count << ", " << sum_delta_element_count/mcomm.get_nprocs() << ") "
                              << (float) max_delta_element_count/min_delta_element_count << " "
                              << std::endl;
                else if (min_full_element_count != 0 && min_delta_element_count == 0)
                    std::cout << "[LOAD BALANCING] Full (" << max_full_element_count << ", " << min_full_element_count << ", " << sum_full_element_count/mcomm.get_nprocs() << ") "
                              << (float) max_full_element_count/min_full_element_count << " "
                              << "Delta (" << max_delta_element_count << ", " << min_delta_element_count << ", " << sum_delta_element_count/mcomm.get_nprocs() << ") "
                              << std::endl;
                else
                    std::cout << "[LOAD BALANCING] Full (" << max_full_element_count << ", " << min_full_element_count << ", " << sum_full_element_count/mcomm.get_nprocs() << ") "
                              << "Delta (" << max_delta_element_count << ", " << min_delta_element_count << ", " << sum_delta_element_count/mcomm.get_nprocs() << ") "
                              << std::endl;
            }
        }

        return;
    }



    void initial_load_balance(float r_factor)
    {
        for (std::vector<relation*>::iterator it = relation_manager.begin() ; it != relation_manager.end(); ++it)
        {
            relation* current_relation = *it;
            if (load_balance_merge_full_and_delta(current_relation, r_factor) == false)
                load_balance_split_full_and_delta(current_relation, r_factor, rc);
        }

        return;
    }
#endif


    void execute()
    {
        u32 local_join_count = 0;
#if 1
        double running_intra_bucket_comm_time = 0;
        double running_local_join_time = 0;
        double running_all_to_all_time = 0;
        double running_insert_in_full_time = 0;
        double running_insert_in_newt_time = 0;
        double running_verify_time = 0;
        double running_lb = 0;

        double intra_bucket_start = 0, intra_bucket_end = 0;
        double local_join_start = 0, local_join_end = 0;
        double all_to_all_start = 0, all_to_all_end = 0;
        double insert_full_start = 0, insert_full_end = 0;
        double insert_newt_start = 0, insert_newt_end = 0;
        double verify_start = 0, verify_end = 0;

        double iteration_time = 0;
        double running_time = 0;
        int outer_loop = 0;
        int inner_loop = 0;
        int iteration = 1;
#endif
        u32 RA_count = RA_list.size();
        int *offset = new int[RA_count];
        for (u32 i =0; i < RA_count; i++)
            offset[i] = 0;

        //bool threshold_reached = false;
        int rank = mcomm.get_rank();


        double start_time = MPI_Wtime();



        //if (rank == 0)
        //    std::cout <<  "Threshold " << threshold << " RF " << refinement_factor << " RI " << refinement_ts << std::endl;

#if 0
        if (refinement_ts != 0)
            if (outer_loop % refinement_ts == 0)
                if (!threshold_reached)
                    load_balance(refinement_factor, refinement_chooser);
#endif

        print_full();
        while (true)
        {
            intra_bucket_start = MPI_Wtime();
            //if (!threshold_reached) intra_bucket_comm();
            intra_bucket_comm();
            intra_bucket_end = MPI_Wtime();
            running_intra_bucket_comm_time = running_intra_bucket_comm_time + (intra_bucket_end - intra_bucket_start);


            local_join_start = MPI_Wtime();
            //threshold_reached = local_compute(&local_join_count, offset);
            local_compute(&local_join_count, offset);
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
            //local_insert_in_full(threshold_reached);
            local_insert_in_full();
            insert_full_end = MPI_Wtime();
            running_insert_in_full_time = running_insert_in_full_time + (insert_full_end - insert_full_start);
#if 1
            /*
            if (threshold_reached == true)
            {
                inner_loop = 0;
                outer_loop++;
            }
            */

            //if (rank == 0)
            //    std::cout << "RUNNING ITERATION [" << iteration <<"] " << running_intra_bucket_comm_count << " " << running_local_join_count << " " << running_all_to_all_count << " " << running_insert_in_full_count << " " << running_insert_in_delta << std::endl;

            verify_start = MPI_Wtime();
            //if (check_for_fixed_point(threshold_reached) == true)
            if (check_for_fixed_point() == true)
            {
                verify_end = MPI_Wtime();
                running_verify_time = running_verify_time + (verify_end - verify_start);
                iteration_time = (verify_end - verify_start) + (insert_newt_end - insert_newt_start) + (insert_full_end - insert_full_start) + (all_to_all_end - all_to_all_start) + (local_join_end - local_join_start) + (intra_bucket_end - intra_bucket_start);
                running_time = running_time + iteration_time;

                if (rank == 0)
                    std::cout << "T CC OL " << outer_loop << " IL " << inner_loop << " [" << iteration << "] " << running_time << " " << iteration_time
                              << " intra_bucket " <<  (intra_bucket_end - intra_bucket_start)
                              << " LJ " <<  (local_join_end - local_join_start)
                              << " All to All " <<  (all_to_all_end - all_to_all_start)
                              << " Insert Full " <<  (insert_full_end - insert_full_start)
                              << " Insert newt " <<  (insert_newt_end - insert_newt_start)
                              << " Verify " <<  (verify_end - verify_start)
                              << std::endl;
                break;
            }

            /*
            parallel_RA* current_r = RA_list[0];
            if (current_r->get_RA_type() == COPY)
            {
                //std::cout << "COPY ONLY" << std::endl;
                break;
            }
            */


            verify_end = MPI_Wtime();
            running_verify_time = running_verify_time + (verify_end - verify_start);

            iteration_time = (verify_end - verify_start) + (insert_newt_end - insert_newt_start) + (insert_full_end - insert_full_start) + (all_to_all_end - all_to_all_start) + (local_join_end - local_join_start) + (intra_bucket_end - intra_bucket_start);
            running_time = running_time + iteration_time;

            if (rank == 0)
                std::cout << "F CC OL " << outer_loop << " IL " << inner_loop << " [" << iteration << "] " << running_time << " " << iteration_time
                          << " intra_bucket " <<  (intra_bucket_end - intra_bucket_start)
                          << " LJ " <<  (local_join_end - local_join_start)
                          << " All to All " <<  (all_to_all_end - all_to_all_start)
                          << " Insert Full " <<  (insert_full_end - insert_full_start)
                          << " Insert newt " <<  (insert_newt_end - insert_newt_start)
                          << " Verify " <<  (verify_end - verify_start)
                          << std::endl;


            inner_loop++;
            iteration++;

            /*
            if (current_r->get_iteration_count() != -1)
            {
                current_r->decrement_iteration_count();
                if (current_r->get_iteration_count() == 0)
                    break;
            }
            */
#endif
        }

#if 1
        double end_time = MPI_Wtime();
        delete[] offset;

        double total_time = end_time - start_time;
        double max_time = 0;
        MPI_Allreduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        if (total_time == max_time)
            std::cout << "Rank " << rank
                      << " Total Time: [" << (end_time - start_time)
                      << " " << running_time
                      << " "
                      << (running_intra_bucket_comm_time + running_local_join_time + running_all_to_all_time + running_insert_in_newt_time + running_insert_in_full_time + running_lb + running_verify_time)
                      << "] intra_bucket " << running_intra_bucket_comm_time
                      << " LJ " << running_local_join_time
                      << " A2A " << running_all_to_all_time
                      << " Insert in new " << running_insert_in_newt_time
                      << " Insert in full " << running_insert_in_full_time
                      << " LB " << running_lb
                      << " FPC " << running_verify_time << std::endl;
        print_full();
#endif
    }

};

#endif
