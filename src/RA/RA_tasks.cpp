#include "../parallel_RA_inc.h"


void RAM::set_comm(mpi_comm& mcomm)
{
    this->mcomm = mcomm;
    for (parallel_RA* ra : RA_list)
        ra->set_comm(mcomm);
    for (relation* r : relation_manager)
        r->set_mcomm(mcomm);
}


// Print relation size
// this function is hacky right now
void print_full(int tid, int loop_count, std::vector<relation*> relation_manager, mpi_comm mcomm)
{
    u64 gfull = 0, tfull = 0;
    u32 counter = 0;
    for (std::vector<relation*>::iterator it = relation_manager.begin() ; it != relation_manager.end(); ++it)
    {
        relation* rel = *it;
        u64 send = rel->get_full_element_count();
        u64 recv;
        MPI_Allreduce(&send, &recv, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mcomm.get_local_comm());

        if (counter == 0)
            gfull = recv;
        else
            tfull = recv;
        counter++;
#if 0
        if (rel->get_filename() != NULL)
        {
            std::vector<u64> prefix = {};
            vector_buffer vb = vector_buffer_create_empty();
            google_relation* full = rel->get_full();

            full[mcomm.get_local_rank()].as_vector_buffer(&vb, prefix, 0);

            if (mcomm.get_local_rank() == 0)
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

            MPI_Barrier(mcomm.get_local_comm());

            int int_send = send;
            int scan;
            MPI_Scan(&int_send, &scan, 1, MPI_INT, MPI_SUM, mcomm.get_local_comm());
            MPI_File fp;
            MPI_Status status;
            MPI_Offset offset = scan;
            char fname[PATH_MAX];
            sprintf(fname, "%s%s", rel->get_filename(), "/data.raw");
            std::cout << "fname " << fname << "Offset " << offset << " send " << int_send << std::endl;
            MPI_File_open(mcomm.get_local_comm(), fname, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
            MPI_File_write_at_all(fp, offset, vb.buffer, int_send, MPI_UNSIGNED_LONG_LONG, &status);
            MPI_File_close(&fp);
        }
#endif
    }

    if (mcomm.get_local_rank() == 0)
        std::cout << "[Task id " << tid << " ] [Iteration count " << loop_count << " [ " << mcomm.get_rank() << " " << mcomm.get_local_nprocs() << "] G " << gfull << " T " << tfull << std::endl;

    return;
}


// check for fixed point
void RAM::check_for_fixed_point(std::vector<u64>& history)
{
    int delta_sum = 0;
    int full_sum = 0;
    for (relation* rel: relation_manager)
    {
        int sum1 = 0;
        int delta_element_count = rel->get_delta_element_count();
        MPI_Allreduce(&delta_element_count, &sum1, 1, MPI_INT, MPI_SUM, mcomm.get_local_comm());
        delta_sum = delta_sum + sum1;

        int sum2 = 0;
        int full_element_count = rel->get_full_element_count();
        MPI_Allreduce(&full_element_count, &sum2, 1, MPI_INT, MPI_SUM, mcomm.get_local_comm());
        full_sum = full_sum + sum2;
    }
    history.push_back(delta_sum);
    history.push_back(full_sum);
}


// intra_bucket Comm
u64 RAM::intra_bucket_comm_execute()
{
    u64 total_data_moved = 0;
    u32 counter = 0;
    u32 RA_count = RA_list.size();

    intra_bucket_buf_output_size = new u64[RA_count];
    intra_bucket_buf_output = new u64*[RA_count];

    for (std::vector<parallel_RA*>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
    {
        // No intra-bucket comm required for copy
        if ((*it)->get_RA_type() == COPY)
            continue;


        // No intra-bucket comm required for acopy
        else if ((*it)->get_RA_type() == ACOPY)
            continue;

        // Intra-bucket comm for joins
        else if ((*it)->get_RA_type() == JOIN)
        {
            parallel_join* current_ra = (parallel_join*) *it;


            // Join between delta and delta
            if (current_ra->get_join_input0_graph_type() == DELTA && current_ra->get_join_input1_graph_type() == DELTA)
            {
                // We are picking the second delta as the source (ROOM FOR OPTIMIZATION HERE)
                relation* input = current_ra->get_join_input1();
                relation* output = current_ra->get_join_input0();

                intra_bucket_comm(get_bucket_count(),
                                              input->get_delta(),
                                              input->get_distinct_sub_bucket_rank_count(), input->get_distinct_sub_bucket_rank(), input->get_bucket_map(),
                                              output->get_distinct_sub_bucket_rank_count(), output->get_distinct_sub_bucket_rank(), output->get_bucket_map(),
                                              &intra_bucket_buf_output_size[counter], &intra_bucket_buf_output[counter],
                                              mcomm.get_local_comm());
                total_data_moved = total_data_moved + intra_bucket_buf_output_size[counter];
            }

            // Join between delta and full
            else if (current_ra->get_join_input0_graph_type() == DELTA && current_ra->get_join_input1_graph_type() == FULL)
            {
                relation* input = current_ra->get_join_input0();
                relation* output = current_ra->get_join_input1();

                intra_bucket_comm(get_bucket_count(),
                                              input->get_delta(),
                                              input->get_distinct_sub_bucket_rank_count(), input->get_distinct_sub_bucket_rank(), input->get_bucket_map(),
                                              output->get_distinct_sub_bucket_rank_count(), output->get_distinct_sub_bucket_rank(), output->get_bucket_map(),
                                              &intra_bucket_buf_output_size[counter], &intra_bucket_buf_output[counter],
                                              mcomm.get_local_comm());
                total_data_moved = total_data_moved + intra_bucket_buf_output_size[counter];

            }

            // Join between full and delta
            else if (current_ra->get_join_input0_graph_type() == FULL && current_ra->get_join_input1_graph_type() == DELTA)
            {
                relation* input = current_ra->get_join_input1();
                relation* output = current_ra->get_join_input0();

                intra_bucket_comm(get_bucket_count(),
                                              input->get_delta(),
                                              input->get_distinct_sub_bucket_rank_count(), input->get_distinct_sub_bucket_rank(), input->get_bucket_map(),
                                              output->get_distinct_sub_bucket_rank_count(), output->get_distinct_sub_bucket_rank(), output->get_bucket_map(),
                                              &intra_bucket_buf_output_size[counter], &intra_bucket_buf_output[counter],
                                              mcomm.get_local_comm());
                total_data_moved = total_data_moved + intra_bucket_buf_output_size[counter];
            }

            // Join between full and full
            else if (current_ra->get_join_input0_graph_type() == FULL && current_ra->get_join_input1_graph_type() == FULL)
            {
                relation* input = current_ra->get_join_input1();
                relation* output = current_ra->get_join_input0();

                intra_bucket_comm(get_bucket_count(),
                                              input->get_full(),
                                              input->get_distinct_sub_bucket_rank_count(), input->get_distinct_sub_bucket_rank(), input->get_bucket_map(),
                                              output->get_distinct_sub_bucket_rank_count(), output->get_distinct_sub_bucket_rank(), output->get_bucket_map(),
                                              &intra_bucket_buf_output_size[counter], &intra_bucket_buf_output[counter],
                                              mcomm.get_local_comm());
                total_data_moved = total_data_moved + intra_bucket_buf_output_size[counter];
            }
        }
        counter++;
    }

    return total_data_moved;
}


u32 RAM::allocate_compute_buffers()
{
    u32 allocated_memory_size = 0;
    compute_buffer.ra_count = RA_list.size();
    compute_buffer.nprocs = get_bucket_count();
    compute_buffer.ra_size = new int[compute_buffer.ra_count];

    allocated_memory_size = allocated_memory_size + compute_buffer.ra_count * sizeof(int);

    compute_buffer.local_compute_output = new vector_buffer*[compute_buffer.ra_count];
    compute_buffer.local_compute_output_size = new int*[compute_buffer.ra_count];

    allocated_memory_size = allocated_memory_size + compute_buffer.ra_count * sizeof(vector_buffer*);
    allocated_memory_size = allocated_memory_size + compute_buffer.ra_count * sizeof(int*);

    compute_buffer.cumulative_tuple_process_map = new int[compute_buffer.nprocs];
    memset(compute_buffer.cumulative_tuple_process_map, 0, compute_buffer.nprocs * sizeof(int));
    allocated_memory_size = allocated_memory_size + compute_buffer.nprocs * sizeof(int);

    for (int i = 0; i < compute_buffer.ra_count; i++)
    {
        compute_buffer.local_compute_output[i] = new vector_buffer[compute_buffer.nprocs];
        allocated_memory_size = allocated_memory_size + compute_buffer.nprocs * sizeof(vector_buffer);

        for (int j = 0; j < compute_buffer.nprocs; j++)
            compute_buffer.local_compute_output[i][j].vector_buffer_create_empty();

        compute_buffer.local_compute_output_size[i] = new int[compute_buffer.nprocs];
        memset(compute_buffer.local_compute_output_size[i], 0, compute_buffer.nprocs * sizeof(int));
        allocated_memory_size = allocated_memory_size + compute_buffer.nprocs * sizeof(int);
    }

    return allocated_memory_size;
}


void RAM::free_compute_buffers()
{
    for (int i = 0; i < compute_buffer.ra_count; i++)
    {
        delete[] compute_buffer.local_compute_output[i];
        delete[] compute_buffer.local_compute_output_size[i];
    }
    delete[] compute_buffer.ra_size;
    delete[] compute_buffer.local_compute_output;
    delete[] compute_buffer.local_compute_output_size;
    delete[] compute_buffer.cumulative_tuple_process_map;

}


u32 RAM::local_compute()
{
    u32 join_tuples = 0;
    u32 join_tuples_duplicates = 0;
    u32 total_join_tuples = 0;
    //u32 RA_count = RA_list.size();
    //u32 nprocs = (u32)mcomm.get_local_nprocs();

    //bool join_completed = true;


    u32 counter = 0;
    //u32 current_join_duplicates = 0;

    for (std::vector<parallel_RA*>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
    {
        if ((*it)->get_RA_type() == COPY)
        {
            parallel_copy* current_ra = (parallel_copy*) *it;


            std::vector<int> reorder_map_array;
            current_ra->get_copy_rename_index(&reorder_map_array);

            if (current_ra->get_copy_input0_graph_type() == DELTA)
            {
                relation* output_relation = current_ra->get_copy_output();
                relation* input_relation = current_ra->get_copy_input();

                current_ra->local_copy(get_bucket_count(),
                                       input_relation->get_delta(), input_relation->get_bucket_map(),
                                       output_relation,
                                       reorder_map_array,
                                       input_relation->get_arity(),
                                       input_relation->get_join_column_count(),
                                       compute_buffer, counter);
            }
            if (current_ra->get_copy_input0_graph_type() == FULL)
            {
                relation* output_relation = current_ra->get_copy_output();
                relation* input_relation = current_ra->get_copy_input();

                current_ra->local_copy(get_bucket_count(),
                                       input_relation->get_full(), input_relation->get_bucket_map(),
                                       output_relation,
                                       reorder_map_array,
                                       input_relation->get_arity(),
                                       input_relation->get_join_column_count(),
                                       compute_buffer, counter);
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


            std::vector<int> reorder_map_array;
            current_ra->get_join_projection_index(&reorder_map_array);

            if (current_ra->get_join_input0_graph_type() == DELTA && current_ra->get_join_input1_graph_type() == DELTA)
            {
                relation* delta0 = current_ra->get_join_input0();
                relation* delta1 = current_ra->get_join_input1();

                //join_completed = join_completed &
                current_ra->local_join(get_bucket_count(),
                                                                 intra_bucket_buf_output_size[counter], delta1->get_arity(), intra_bucket_buf_output[counter], 0,
                                                                 delta0->get_delta(), delta0->get_delta_element_count(), delta0->get_arity(),
                                                                 reorder_map_array,
                                                                 output_relation,
                                                                 compute_buffer,
                                                                 counter,
                                                                 join_column_count,
                                                                 &join_tuples_duplicates,
                                                                 &join_tuples);
                total_join_tuples = total_join_tuples + join_tuples;
            }
            else if (current_ra->get_join_input0_graph_type() == DELTA && current_ra->get_join_input1_graph_type() == FULL)
            {
                relation* full = current_ra->get_join_input1();
                relation* delta = current_ra->get_join_input0();


                //join_completed = join_completed &
                current_ra->local_join(get_bucket_count(),
                                                                 intra_bucket_buf_output_size[counter], delta->get_arity(), intra_bucket_buf_output[counter], 0,
                                                                 full->get_full(), full->get_full_element_count(), full->get_arity(),
                                                                 reorder_map_array,
                                                                 output_relation,
                                                                 compute_buffer,
                                                                 counter,
                                                                 join_column_count,
                                                                 &join_tuples_duplicates,
                                                                 &join_tuples);
                total_join_tuples = total_join_tuples + join_tuples;
            }
            else if (current_ra->get_join_input0_graph_type() == FULL && current_ra->get_join_input1_graph_type() == DELTA)
            {
                relation* full = current_ra->get_join_input0();
                relation* delta = current_ra->get_join_input1();

                //join_completed = join_completed &
                current_ra->local_join(get_bucket_count(),
                                                                 intra_bucket_buf_output_size[counter], delta->get_arity(), intra_bucket_buf_output[counter], 1,
                                                                 full->get_full(), full->get_full_element_count(), full->get_arity(),
                                                                 reorder_map_array,
                                                                 output_relation,
                                                                 compute_buffer,
                                                                 counter,
                                                                 join_column_count,
                                                                 &join_tuples_duplicates,
                                                                 &join_tuples);
                total_join_tuples = total_join_tuples + join_tuples;
            }
        }
        counter++;
    }

    counter = 0;
    for (std::vector<parallel_RA*>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
    {
        parallel_RA* current_ra = *it;
        if (current_ra->get_RA_type() == JOIN)
            delete[] intra_bucket_buf_output[counter];
    }

    delete[] intra_bucket_buf_output_size;
    delete[] intra_bucket_buf_output;


    return join_tuples;
}



int RAM::all_to_all()
{
    int counter = 0;
    u64 total_all_to_all = 0;
    u32 RA_count = compute_buffer.ra_count;
    int nprocs = compute_buffer.nprocs;

    cumulative_all_to_all_recv_process_size_array = new int[RA_count * nprocs];

    int *flat_process_size = new int[nprocs * RA_count];
    for (int i = 0; i < nprocs; i++)
        for (u32 j = 0; j < RA_count; j++)
            flat_process_size[counter++] = compute_buffer.local_compute_output_size[j][i];

    /* This step prepares for actual data transfer */
    /* Every process sends to every other process the amount of data it is going to send */
    //int recv_process_size_array[nprocs * RA_count];
    memset(cumulative_all_to_all_recv_process_size_array, 0, RA_count * nprocs * sizeof(int));
    MPI_Alltoall(flat_process_size, RA_count, MPI_INT, cumulative_all_to_all_recv_process_size_array, RA_count, MPI_INT, mcomm.get_local_comm());

    delete[] flat_process_size;


    // PREPARING THE SEND PREFIX
    int cumulative_prefix_sum_process_size[nprocs];
    memset(cumulative_prefix_sum_process_size, 0, nprocs * sizeof(int));

    for (int i = 1; i < nprocs; i++)
        cumulative_prefix_sum_process_size[i] = cumulative_prefix_sum_process_size[i - 1] + compute_buffer.cumulative_tuple_process_map[i - 1];

    int process_data_buffer_size = cumulative_prefix_sum_process_size[nprocs - 1] + compute_buffer.cumulative_tuple_process_map[nprocs - 1];


    // PREPARING THE SEND DATA
    u64* process_data = 0;
    process_data = new u64[process_data_buffer_size];
    memset(process_data, 0, process_data_buffer_size * sizeof(u64));

    u32 boffset = 0;
    for(int i = 0; i < nprocs; i++)
    {
        for (u32 r = 0; r < RA_count; r++)
        {
            memcpy(process_data + boffset, compute_buffer.local_compute_output[r][i].buffer, compute_buffer.local_compute_output[r][i].size);
            boffset = boffset + (compute_buffer.local_compute_output[r][i].size)/8;
            compute_buffer.local_compute_output[r][i].vector_buffer_free();
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
    MPI_Alltoallv(process_data, compute_buffer.cumulative_tuple_process_map, cumulative_prefix_sum_process_size, MPI_UNSIGNED_LONG_LONG, cumulative_all_to_all_buffer, cumulative_recv_process_size_array, cumulative_prefix_sum_recv_process_size_buffer, MPI_UNSIGNED_LONG_LONG, mcomm.get_local_comm());

    delete[] process_data;

#if 0
    for (u32 i = 0; i < RA_count; i++)
    {
        delete[] local_compute_output[i];
        delete[] local_compute_output_size[i];
    }

    delete[] local_compute_output;
    delete[] local_compute_output_size;
    delete[] cumulative_tuple_process_map;
#endif


    return total_all_to_all;
}


void RAM::local_insert_in_newt()
{
    u32 successful_insert = 0;
    int nprocs = mcomm.get_local_nprocs();
    int RA_count = RA_list.size();
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
        for (u32 x = starting; x < starting + elements_to_read; x=x+arity)
        {
            u64 t[arity];
            for (u32 a = 0; a < arity; a++)
                t[a] = cumulative_all_to_all_buffer[x + a];

            if (output->find_in_full(t) == false)
                if (output->insert_in_newt(t) == true)
                    successful_insert++;
        }
        starting = starting + elements_to_read;
    }

    delete[] cumulative_all_to_all_recv_process_size_array;
    delete[] cumulative_all_to_all_buffer;
}


void RAM::local_insert_in_full()
{
    for (std::vector<relation*>::iterator it = relation_manager.begin() ; it != relation_manager.end(); ++it)
    {
        relation* current_r = *it;
        current_r->insert_delta_in_full();
        current_r->local_insert_in_delta();
    }
    return;
}



void RAM::load_balance()
{
    for (std::vector<relation*>::iterator it = relation_manager.begin() ; it != relation_manager.end(); ++it)
    {
        relation* current_relation = *it;
        if (current_relation->load_balance_merge_full_and_delta(refinement_factor) == false)
            current_relation->load_balance_split_full_and_delta(refinement_factor);

        u64 max_full_element_count = 0;
        u64 min_full_element_count = 0;
        u64 sum_full_element_count = 0;
        u64 full_element_count = 0;

        u64 max_delta_element_count = 0;
        u64 min_delta_element_count = 0;
        u64 sum_delta_element_count = 0;
        u64 delta_element_count = 0;

        u32* bucket_map = current_relation->get_bucket_map();
        u32* sub_bucket_count = current_relation->get_sub_bucket_per_bucket_count();

        u32** full_sub_bucket_size = current_relation->get_full_sub_bucket_element_count();
        u32** delta_sub_bucket_size = current_relation->get_delta_sub_bucket_element_count();

        for (u32 i = 0; i < get_bucket_count(); i++)
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

        MPI_Allreduce(&full_element_count, &max_full_element_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, mcomm.get_local_comm());
        MPI_Allreduce(&full_element_count, &min_full_element_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, mcomm.get_local_comm());
        MPI_Allreduce(&full_element_count, &sum_full_element_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mcomm.get_local_comm());

        MPI_Allreduce(&delta_element_count, &max_delta_element_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, mcomm.get_local_comm());
        MPI_Allreduce(&delta_element_count, &min_delta_element_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, mcomm.get_local_comm());
        MPI_Allreduce(&delta_element_count, &sum_delta_element_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mcomm.get_local_comm());

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



void RAM::print_all_relation()
{
    u32 rel_count=0;
    for (std::vector<relation*>::iterator it = relation_manager.begin() ; it != relation_manager.end(); ++it)
    {
        std::cout << "Relation " << rel_count++ << std::endl;

        relation* current_r = *it;
        current_r->print();

        std::cout << std::endl;
    }
    return;
}



void RAM::enable_logging()
{
    logging = true;
    return;
}



void RAM::execute(int batch_size, std::vector<u64>& history, int task_id)
{
    if (mcomm.get_rank() == 0)
        std::cout << std::endl <<  "-------------------------------------------------------------------------------------------------------------------------------" << std::endl;

    int inner_loop = 0;
    double running_time = 0;
    while (batch_size != 0)
    {   
        double intra_start = MPI_Wtime();
        intra_bucket_comm_execute();
        double intra_end = MPI_Wtime();

        allocate_compute_buffers();

        double compute_start = MPI_Wtime();
        local_compute();
        double compute_end = MPI_Wtime();

        double all_to_all_start = MPI_Wtime();
        all_to_all();
        double all_to_all_end = MPI_Wtime();

        free_compute_buffers();

        double insert_in_newt_start = MPI_Wtime();
        local_insert_in_newt();
        double insert_in_newt_end = MPI_Wtime();

        double insert_in_full_start = MPI_Wtime();
        local_insert_in_full();
        double insert_in_full_end = MPI_Wtime();

        running_time = running_time + (intra_end - intra_start) + (compute_end - compute_start) + (all_to_all_end - all_to_all_start) + (insert_in_newt_end - insert_in_newt_start) + (insert_in_full_end - insert_in_full_start);

#if 0
        if (mcomm.get_rank() == 0)
            std::cout << "INNER [" << inner_loop << "] "
                      << " Intra " << (intra_end - intra_start)
                      << " compute " << (compute_end - compute_start)
                      << " All to all " << (all_to_all_end - all_to_all_start)
                      << " Insert newt " << (insert_in_newt_end - insert_in_newt_start)
                      << " Insert full " << (insert_in_full_end - insert_in_full_start)
                      << " Total " << (intra_end - intra_start) + (compute_end - compute_start) + (all_to_all_end - all_to_all_start) + (insert_in_newt_end - insert_in_newt_start) + (insert_in_full_end - insert_in_full_start)
                      << " [ "
                      << running_time
                      << " ]" << std::endl;
#endif

        batch_size--;
        inner_loop++;

        if (iteration_count == 1)
            break;
    }

    check_for_fixed_point(history);

    if (enable_dump_io == true)
        print_full(task_id, inner_loop, relation_manager, mcomm);
}



void RAM::execute_time(double batch_time, std::vector<u64>& history, int task_id)
{
    double running_total = 0;
    double running_time = 0;
    int inner_loop = 0;

    if (mcomm.get_rank() == 0)
        std::cout << std::endl <<  "-------------------------------------------------------------------------------------------------------------------------------" << std::endl;

    while (running_time <= batch_time)
    {
        double elapsed_start = MPI_Wtime();

        double intra_start = MPI_Wtime();
        intra_bucket_comm_execute();
        double intra_end = MPI_Wtime();

        allocate_compute_buffers();

        double compute_start = MPI_Wtime();
        local_compute();
        double compute_end = MPI_Wtime();

        double all_to_all_start = MPI_Wtime();
        all_to_all();
        double all_to_all_end = MPI_Wtime();

        free_compute_buffers();

        double insert_in_newt_start = MPI_Wtime();
        local_insert_in_newt();
        double insert_in_newt_end = MPI_Wtime();

        double insert_in_full_start = MPI_Wtime();
        local_insert_in_full();
        double insert_in_full_end = MPI_Wtime();


        running_total = running_total + (intra_end - intra_start) + (compute_end - compute_start) + (all_to_all_end - all_to_all_start) + (insert_in_newt_end - insert_in_newt_start) + (insert_in_full_end - insert_in_full_start);

#if 0
        if (mcomm.get_rank() == 0)
            std::cout << "INNER [" << inner_loop << "] "
                      << " Intra " << (intra_end - intra_start)
                      << " compute " << (compute_end - compute_start)
                      << " All to all " << (all_to_all_end - all_to_all_start)
                      << " Insert newt " << (insert_in_newt_end - insert_in_newt_start)
                      << " Insert full " << (insert_in_full_end - insert_in_full_start)
                      << " Total " << (intra_end - intra_start) + (compute_end - compute_start) + (all_to_all_end - all_to_all_start) + (insert_in_newt_end - insert_in_newt_start) + (insert_in_full_end - insert_in_full_start)
                      << " [ "
                      << running_time
                      << " ]" << std::endl;
#endif
        inner_loop++;

        if (iteration_count == 1)
            break;

        double elapsed_end = MPI_Wtime();
        running_time = running_time + (elapsed_end - elapsed_start);

        double mint;
        MPI_Allreduce(&running_time, &mint, 1, MPI_DOUBLE, MPI_MAX, mcomm.get_local_comm());
        running_time = mint;
    }

    check_for_fixed_point(history);

    if (enable_dump_io == true)
        print_full(task_id, inner_loop, relation_manager, mcomm);
}
