/*
 * scc (tasks)
 * Copyright (c) Sidharth Kumar, et al, see License.md
 */



#include "../parallel_RA_inc.h"

RAM::~RAM()
{
    loop_count_tracker = 0;
    for (std::vector<parallel_RA*>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
        delete (*it);
}

/// Example: RAM* scc13237 = new RAM(true, 1);
/// true: run this scc till fixed point is reached (false: run this scc for only one iteration, for copy and acopy rules)
/// 1: id of scc (used for internal debugging)
RAM::RAM (bool ic, int r_id)
{
    ram_relation_count = 0;
    loop_count_tracker = 0;
    if (ic==false)
        iteration_count=1;
    ram_id = r_id;
    RA_list = {};
}


/// Example
/// rel_edge_2_2: relation that is being added to the SCC
/// false: keep the delta and full the way they are (true: move whatever is in full to delta, once before the start of the fixed point loop)
void RAM::add_relation(relation*& G, bool i_status)
{
    ram_relations[ram_relation_count] = G;
    ram_relation_status[ram_relation_count] = i_status;
    ram_relation_count++;
}



void RAM::set_comm(mpi_comm& mcomm)
{
    this->mcomm = mcomm;

    for (parallel_RA* ra : RA_list)
        ra->set_comm(mcomm);

    for (u32 i=0; i < ram_relation_count; i++)
        ram_relations[i]->set_mcomm(mcomm);
}



void RAM::print_all_relation()
{
    for (u32 i=0; i < ram_relation_count; i++)
        ram_relations[i]->print();
}




void RAM::load_balance()
{
    for (u32 i=0; i < ram_relation_count; i++)
    {
        relation* current_relation = ram_relations[i];
        if (current_relation->load_balance_merge_full_and_delta(refinement_factor) == false)
            current_relation->load_balance_split_full_and_delta(refinement_factor);

        u64 max_full_element_count = 0, min_full_element_count = 0, sum_full_element_count = 0, full_element_count = 0, max_delta_element_count = 0, min_delta_element_count = 0, sum_delta_element_count = 0, delta_element_count = 0;

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
}




u64 RAM::intra_bucket_comm_execute()
{
    u64 total_data_moved = 0;
    u32 counter = 0;
    u32 RA_count = RA_list.size();

    intra_bucket_buf_output_size = new u64[RA_count];
    memset(intra_bucket_buf_output_size, 0, RA_count * sizeof(u64));

    intra_bucket_buf_output = new u64*[RA_count];

    for (std::vector<parallel_RA*>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
    {
        /// No intra-bucket comm required for copy
        if ((*it)->get_RA_type() == COPY)
        {
            counter++;
            continue;
        }

        /// No intra-bucket comm required for copy
        if ((*it)->get_RA_type() == COPY_FILTER)
        {
            counter++;
            continue;
        }

        /// No intra-bucket comm required for acopy
        else if ((*it)->get_RA_type() == ACOPY)
        {
            counter++;
            continue;
        }

        /// Intra-bucket comm for joins
        else if ((*it)->get_RA_type() == JOIN)
        {
            parallel_join* current_ra = (parallel_join*) *it;
            relation* input0 = current_ra->get_join_input0();
            relation* input1 = current_ra->get_join_input1();

            /// Join between delta and delta
            if (current_ra->get_join_input0_graph_type() == DELTA && current_ra->get_join_input1_graph_type() == DELTA)
            {

                intra_bucket_comm(get_bucket_count(),
                                  input0->get_delta(),
                                  input0->get_distinct_sub_bucket_rank_count(), input0->get_distinct_sub_bucket_rank(), input0->get_bucket_map(),
                                  input1->get_distinct_sub_bucket_rank_count(), input1->get_distinct_sub_bucket_rank(), input1->get_bucket_map(),
                                  &intra_bucket_buf_output_size[counter], &intra_bucket_buf_output[counter],
                                  mcomm.get_local_comm());

                total_data_moved = total_data_moved + intra_bucket_buf_output_size[counter];
            }

            /// Join between delta and full
            else if (current_ra->get_join_input0_graph_type() == DELTA && current_ra->get_join_input1_graph_type() == FULL)
            {

                intra_bucket_comm(get_bucket_count(),
                                  input0->get_delta(),
                                  input0->get_distinct_sub_bucket_rank_count(), input0->get_distinct_sub_bucket_rank(), input0->get_bucket_map(),
                                  input1->get_distinct_sub_bucket_rank_count(), input1->get_distinct_sub_bucket_rank(), input1->get_bucket_map(),
                                  &intra_bucket_buf_output_size[counter], &intra_bucket_buf_output[counter],
                                  mcomm.get_local_comm());
                total_data_moved = total_data_moved + intra_bucket_buf_output_size[counter];
            }

            /// Join between full and delta
            else if (current_ra->get_join_input0_graph_type() == FULL && current_ra->get_join_input1_graph_type() == DELTA)
            {
                intra_bucket_comm(get_bucket_count(),
                                  input1->get_delta(),
                                  input1->get_distinct_sub_bucket_rank_count(), input1->get_distinct_sub_bucket_rank(), input1->get_bucket_map(),
                                  input0->get_distinct_sub_bucket_rank_count(), input0->get_distinct_sub_bucket_rank(), input0->get_bucket_map(),
                                  &intra_bucket_buf_output_size[counter], &intra_bucket_buf_output[counter],
                                  mcomm.get_local_comm());
                total_data_moved = total_data_moved + intra_bucket_buf_output_size[counter];
            }

            /// Join between full and full
            else if (current_ra->get_join_input0_graph_type() == FULL && current_ra->get_join_input1_graph_type() == FULL)
            {

                intra_bucket_comm(get_bucket_count(),
                                  input1->get_full(),
                                  input1->get_distinct_sub_bucket_rank_count(), input1->get_distinct_sub_bucket_rank(), input1->get_bucket_map(),
                                  input0->get_distinct_sub_bucket_rank_count(), input0->get_distinct_sub_bucket_rank(), input0->get_bucket_map(),
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
    compute_buffer.local_compute_output_size_total = 0;

    compute_buffer.width = new int[compute_buffer.ra_count];

    allocated_memory_size = allocated_memory_size + compute_buffer.ra_count * sizeof(int);

    compute_buffer.local_compute_output = new vector_buffer*[compute_buffer.ra_count];
    compute_buffer.local_compute_output_size = new int*[compute_buffer.ra_count];
    compute_buffer.local_compute_output_size_flat = new int[compute_buffer.ra_count * compute_buffer.nprocs];
    memset(compute_buffer.local_compute_output_size_flat, 0, compute_buffer.ra_count * compute_buffer.nprocs * sizeof(int));

    allocated_memory_size = allocated_memory_size + compute_buffer.ra_count * sizeof(vector_buffer*);
    allocated_memory_size = allocated_memory_size + compute_buffer.ra_count * sizeof(int*);
    allocated_memory_size = allocated_memory_size + compute_buffer.ra_count * compute_buffer.nprocs * sizeof(int);

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




u32 RAM::local_compute(int* offset)
{
    bool join_completed = true;
    u32 join_tuples = 0;
    u32 join_tuples_duplicates = 0;
    u32 total_join_tuples = 0;
    u32 counter = 0;
    int threshold = 1048576;

    for (std::vector<parallel_RA*>::iterator it = RA_list.begin() ; it != RA_list.end(); ++it)
    {
        if ((*it)->get_RA_type() == COPY)
        {
            parallel_copy* current_ra = (parallel_copy*) *it;

            std::vector<int> reorder_map_array;
            current_ra->get_copy_rename_index(&reorder_map_array);
            relation* output_relation = current_ra->get_copy_output();
            relation* input_relation = current_ra->get_copy_input();

            if (current_ra->get_copy_input0_graph_type() == DELTA)
            {
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
                current_ra->local_copy(get_bucket_count(),
                                       input_relation->get_full(), input_relation->get_bucket_map(),
                                       output_relation,
                                       reorder_map_array,
                                       input_relation->get_arity(),
                                       input_relation->get_join_column_count(),
                                       compute_buffer, counter);
            }
        }

        else if ((*it)->get_RA_type() == COPY_FILTER)
        {
            parallel_copy_filter* current_ra = (parallel_copy_filter*) *it;

            std::vector<int> reorder_map_array;
            current_ra->get_copy_filter_rename_index(&reorder_map_array);
            relation* output_relation = current_ra->get_copy_filter_output();
            relation* input_relation = current_ra->get_copy_filter_input();

            if (current_ra->get_copy_filter_input0_graph_type() == DELTA)
            {
                current_ra->local_copy_filter(get_bucket_count(),
                                              input_relation->get_delta(), input_relation->get_bucket_map(),
                                              output_relation,
                                              reorder_map_array,
                                              input_relation->get_arity(),
                                              input_relation->get_join_column_count(),
                                              compute_buffer, counter);
            }
            if (current_ra->get_copy_filter_input0_graph_type() == FULL)
            {
                current_ra->local_copy_filter(get_bucket_count(),
                                              input_relation->get_full(), input_relation->get_bucket_map(),
                                              output_relation,
                                              reorder_map_array,
                                              input_relation->get_arity(),
                                              input_relation->get_join_column_count(),
                                              compute_buffer, counter);
            }
        }

        else if ((*it)->get_RA_type() == ACOPY)
        {
            parallel_acopy* current_ra = (parallel_acopy*) *it;

            std::vector<int> reorder_map_array;
            current_ra->get_acopy_rename_index(&reorder_map_array);
            relation* output_relation = current_ra->get_acopy_output();
            relation* input_relation = current_ra->get_acopy_input();

            if (current_ra->get_acopy_input0_graph_type() == DELTA)
            {
                current_ra->local_acopy(get_bucket_count(),
                                        input_relation->get_delta(), input_relation->get_bucket_map(),
                                        output_relation,
                                        reorder_map_array,
                                        input_relation->get_arity(),
                                        input_relation->get_join_column_count(),
                                        compute_buffer, counter);
            }
            else if (current_ra->get_acopy_input0_graph_type() == FULL)
            {
                current_ra->local_acopy(get_bucket_count(),
                                        input_relation->get_full(), input_relation->get_bucket_map(),
                                        output_relation,
                                        reorder_map_array,
                                        input_relation->get_arity(),
                                        input_relation->get_join_column_count(),
                                        compute_buffer, counter);
            }
        }

        else if ((*it)->get_RA_type() == JOIN)
        {
            parallel_join* current_ra = (parallel_join*) *it;
            relation* output_relation = current_ra->get_join_output();

            std::vector<int> reorder_map_array;
            current_ra->get_join_projection_index(&reorder_map_array);
            relation* input0 = current_ra->get_join_input0();
            relation* input1 = current_ra->get_join_input1();
            relation* output = current_ra->get_join_output();
            assert(output->get_arity() == reorder_map_array.size());
            int join_column_count = input0->get_join_column_count();

            if (current_ra->get_join_input0_graph_type() == DELTA && current_ra->get_join_input1_graph_type() == DELTA)
            {
                join_completed = join_completed & current_ra->local_join(threshold, &(offset[counter]),
                                                                         LEFT,
                                                                         get_bucket_count(),
                                                                         intra_bucket_buf_output_size[counter], input0->get_arity()+1, intra_bucket_buf_output[counter],
                                                                         input1->get_delta(), input1->get_delta_element_count(), input1->get_arity()+1,
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

                join_completed = join_completed & current_ra->local_join(threshold, &(offset[counter]),
                                                                         LEFT,
                                                                         get_bucket_count(),
                                                                         intra_bucket_buf_output_size[counter], input0->get_arity()+1, intra_bucket_buf_output[counter],
                                                                         input1->get_full(), input1->get_full_element_count(), input1->get_arity()+1,
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

                join_completed = join_completed & current_ra->local_join(threshold, &(offset[counter]),
                                                                         RIGHT,
                                                                         get_bucket_count(),
                                                                         intra_bucket_buf_output_size[counter], input1->get_arity()+1, intra_bucket_buf_output[counter],
                                                                         input0->get_full(), input0->get_full_element_count(), input0->get_arity()+1,
                                                                         reorder_map_array,
                                                                         output_relation,
                                                                         compute_buffer,
                                                                         counter,
                                                                         join_column_count,
                                                                         &join_tuples_duplicates,
                                                                         &join_tuples);
                total_join_tuples = total_join_tuples + join_tuples;
            }
            else if (current_ra->get_join_input0_graph_type() == FULL && current_ra->get_join_input1_graph_type() == FULL)
            {

                join_completed = join_completed & current_ra->local_join(threshold, &(offset[counter]),
                                                                         RIGHT,
                                                                         get_bucket_count(),
                                                                         intra_bucket_buf_output_size[counter], input1->get_arity()+1, intra_bucket_buf_output[counter],
                                                                         input0->get_full(), input0->get_full_element_count(), input0->get_arity()+1,
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

#if 0
    int global_total_join_tuples = 0;
    int global_join_tuples_duplicates = 0;
    MPI_Allreduce(&total_join_tuples, &global_total_join_tuples, 1, MPI_INT, MPI_SUM, mcomm.get_local_comm());
    MPI_Allreduce(&join_tuples_duplicates, &global_join_tuples_duplicates, 1, MPI_INT, MPI_SUM, mcomm.get_local_comm());
    if (mcomm.get_rank() == 0)
        std::cout << "Joins: " << global_total_join_tuples << " Duplicates " << global_join_tuples_duplicates << " ";
#endif

    int global_synchronizer = 0;
    int synchronizer = 0;
    if (join_completed == true)
        synchronizer = 1;

    MPI_Allreduce(&synchronizer, &global_synchronizer, 1, MPI_INT, MPI_BAND, mcomm.get_comm());
    if (global_synchronizer == 1)
    {
        counter = 0;
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
}



void RAM::all_to_all()
{
    comm_compaction_all_to_all(compute_buffer, &cumulative_all_to_all_recv_process_size_array, &cumulative_all_to_all_buffer, mcomm.get_local_comm());
}



void RAM::free_compute_buffers()
{
    for (int i = 0; i < compute_buffer.ra_count; i++)
    {
        delete[] compute_buffer.local_compute_output[i];
        delete[] compute_buffer.local_compute_output_size[i];
    }
    delete[] compute_buffer.local_compute_output_size_flat;
    delete[] compute_buffer.width;
    delete[] compute_buffer.local_compute_output;
    delete[] compute_buffer.local_compute_output_size;
    delete[] compute_buffer.cumulative_tuple_process_map;
}



void RAM::local_insert_in_newt(std::map<u64, u64>& intern_map)
{
    u32 successful_insert = 0, starting = 0;
    int nprocs = mcomm.get_local_nprocs();
    int RA_count = RA_list.size();
    u64 relation_id=0, bucket_id=0, intern_key=0, intern_value=0;

    for (int k = 0; k < RA_count * nprocs; k++)
    {
        successful_insert = 0;
        u32 ra_id = k % RA_count;
        u32 elements_to_read = cumulative_all_to_all_recv_process_size_array[k];
        relation* output;

        if (RA_list[ra_id]->get_RA_type() == COPY)
            output = RA_list[ra_id]->get_copy_output();
        else if (RA_list[ra_id]->get_RA_type() == COPY_FILTER)
            output = RA_list[ra_id]->get_copy_filter_output();
        else if (RA_list[ra_id]->get_RA_type() == JOIN)
            output = RA_list[ra_id]->get_join_output();
        else
            output = RA_list[ra_id]->get_acopy_output();

        if (RA_list[ra_id]->get_RA_type() == COPY || RA_list[ra_id]->get_RA_type() == JOIN || RA_list[ra_id]->get_RA_type() == COPY_FILTER)
        {
            u32 width = output->get_arity();
            u64 tuple[width + 1];

            for (u32 x = starting; x < starting + elements_to_read; x = x + width)
            {
                if (output->find_in_full(cumulative_all_to_all_buffer + x, width) == false &&
                        output->find_in_delta(cumulative_all_to_all_buffer + x, width) == false &&
                        output->find_in_newt(cumulative_all_to_all_buffer + x, width) == false)
                {
                    for (u32 i = 0; i < width; i++)
                        tuple[i] = cumulative_all_to_all_buffer[x+i];

                    relation_id = output->get_intern_tag();
                    relation_id = relation_id<<46;
                    bucket_id = tuple_hash(tuple, output->get_join_column_count()) % get_bucket_count();
                    bucket_id = bucket_id<<28;

                    intern_key = relation_id | bucket_id;

                    std::map<u64 ,u64>::const_iterator it = intern_map.find(intern_key);
                    if( it == intern_map.end() )
                        intern_value=0;
                    else
                        intern_value = it->second + 1;

                    intern_map[intern_key] = intern_value;
                    tuple[width] = intern_key | intern_value;    /// Intern here

                    if (output->insert_in_newt(tuple) == true)
                        successful_insert++;
                }
            }
        }
        else
        {
            u32 width = output->get_arity() + 1;
            u64 tuple[width];
            successful_insert = 0;
            for (u32 x = starting; x < starting + elements_to_read; x = x + width)
            {
                if (output->find_in_full(cumulative_all_to_all_buffer + x, width) == false && output->find_in_delta(cumulative_all_to_all_buffer + x, width) == false)
                {
                    for (u32 i = 0; i < width; i++)
                        tuple[i] = cumulative_all_to_all_buffer[x+i];

                    if (output->insert_in_newt(tuple) == true)
                        successful_insert++;
                }
            }
        }
        starting = starting + elements_to_read;
        //std::cout << output->get_debug_id() << " successful insert " << successful_insert << std::endl;
    }

    delete[] cumulative_all_to_all_recv_process_size_array;
    delete[] cumulative_all_to_all_buffer;
}



void RAM::local_insert_in_full()
{
    for (u32 i=0; i < ram_relation_count; i++)
        //for (std::map<relation*, bool>::iterator it = ram_relations.begin() ; it != ram_relations.end(); ++it)
    {
        //relation* current_r = it->first;
        relation* current_r = ram_relations[i];
        current_r->insert_delta_in_full();
        current_r->local_insert_in_delta();

        //if (current_r->get_debug_id() == 11)
        //    current_r->print();
    }
    return;
}



void RAM::insert_delta_in_full()
{
    for (u32 i=0; i < ram_relation_count; i++)
        ram_relations[i]->insert_delta_in_full();

    return;
}



void RAM::check_for_fixed_point(std::vector<u32>& history)
{
    int local_delta_sum = 0, local_full_sum = 0, global_delta_sum = 0, global_full_sum = 0;
    for (u32 i=0; i < ram_relation_count; i++)
    {
        local_delta_sum = local_delta_sum + ram_relations[i]->get_delta_element_count();
        local_full_sum = local_full_sum + ram_relations[i]->get_full_element_count();
    }
    MPI_Allreduce(&local_delta_sum, &global_delta_sum, 1, MPI_INT, MPI_SUM, mcomm.get_local_comm());
    MPI_Allreduce(&local_full_sum, &global_full_sum, 1, MPI_INT, MPI_SUM, mcomm.get_local_comm());

    history.push_back(global_delta_sum);
    history.push_back(global_full_sum);
}



void RAM::io_all_relation(int status)
{

    char scc_name[1024];

    if (status == 1)
    {
        sprintf(scc_name, "output/scc-%d-iteration_%d", ram_id, loop_count_tracker);
        mkdir(scc_name, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        for (u32 i = 0 ; i < ram_relation_count; i++)
            ram_relations[i]->serial_IO(scc_name);
    }
    else if (status == 0)
    {
        sprintf(scc_name, "output/scc-%d-initial-facts", ram_id);
        mkdir(scc_name, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        for (u32 i = 0 ; i < ram_relation_count; i++)
            ram_relations[i]->serial_IO(scc_name);
    }
    else if (status == 2)
    {
        sprintf(scc_name, "output/scc-%d-output-facts", ram_id);
        mkdir(scc_name, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        for (u32 i = 0 ; i < ram_relation_count; i++)
            ram_relations[i]->serial_IO(scc_name);
    }
}



void RAM::execute_in_batches(int batch_size, std::vector<u32>& history, std::map<u64, u64>& intern_map, double *running_time, double *running_intra_bucket_comm, double *running_buffer_allocate, double *running_local_compute, double *running_all_to_all, double *running_buffer_free, double *running_insert_newt, double *running_insert_in_full)
{
    int inner_loop = 0;
    u32 RA_count = RA_list.size();

    int *offset = new int[RA_count];
    for (u32 i =0; i < RA_count; i++)
        offset[i] = 0;

    while (batch_size != 0)
    {
#if DEBUG_OUTPUT
        if (mcomm.get_rank() == 0)
            std::cout << "--------------FIXED POINT ITERATION " << loop_count_tracker << "--------------" << std::endl;
#endif


        double intra_start = MPI_Wtime();
        intra_bucket_comm_execute();
        double intra_end = MPI_Wtime();
        *running_intra_bucket_comm = *running_intra_bucket_comm + (intra_end - intra_start);



        bool local_join_status = false;
        while (local_join_status == false)
        {

            double allocate_buffers_start = MPI_Wtime();
            allocate_compute_buffers();
            double allocate_buffers_end = MPI_Wtime();
            *running_buffer_allocate = *running_buffer_allocate + (allocate_buffers_end - allocate_buffers_start);


            double compute_start = MPI_Wtime();
            local_join_status = local_compute(offset);
            double compute_end = MPI_Wtime();
            *running_local_compute = *running_local_compute + (compute_end - compute_start);


            double all_to_all_start = MPI_Wtime();
            all_to_all();
            double all_to_all_end = MPI_Wtime();
            *running_all_to_all = *running_all_to_all + (all_to_all_end - all_to_all_start);


            double free_buffers_start = MPI_Wtime();
            free_compute_buffers();
            double free_buffers_end = MPI_Wtime();
            *running_buffer_free = *running_buffer_free + (free_buffers_end - free_buffers_start);


            double insert_in_newt_start = MPI_Wtime();
            local_insert_in_newt(intern_map);
            double insert_in_newt_end = MPI_Wtime();
            *running_insert_newt = *running_insert_newt + (insert_in_newt_end - insert_in_newt_start);


#if DEBUG_OUTPUT
            if (mcomm.get_rank() == 0)
            {
                std::cout << "Current time INNER LOOP [" << loop_count_tracker << " " << inner_loop << "] "
                          << " Buf cre " << (allocate_buffers_end - allocate_buffers_start)
                          << " comp " << (compute_end - compute_start)
                          << " A2A " << (all_to_all_end - all_to_all_start)
                          << " Buf free " << (free_buffers_end - free_buffers_start)
                          << " newt " << (insert_in_newt_end - insert_in_newt_start)
                          << std::endl;

                std::cout << "Running time INNER LOOP [" << loop_count_tracker << " " << inner_loop << "] "
                          << " Buf cre " << *running_buffer_allocate
                          << " comp " << *running_local_compute
                          << " A2A " << *running_all_to_all
                          << " Buf free " << *running_buffer_free
                          << " newt " << *running_insert_newt
                          << std::endl;
            }
#endif
            inner_loop++;
        }

        double insert_in_full_start = MPI_Wtime();
        local_insert_in_full();
        double insert_in_full_end = MPI_Wtime();

        *running_insert_in_full = *running_insert_in_full + (insert_in_full_end - insert_in_full_start);
        *running_time = *running_time + (insert_in_full_end - intra_start);

#if DEBUG_OUTPUT
        if (mcomm.get_rank() == 0)
        {
            std::cout << "Current time OUTER LOOP [" << loop_count_tracker << " ] "
                      << " Intra " << (intra_end - intra_start)
                      << " full " << (insert_in_full_end - insert_in_full_start)
                      << " Total " << (insert_in_full_end - intra_start)
                      << " [ "
                      << *running_time
                      << " ]" << std::endl;

            std::cout << "Running time OUTER LOOP [" << loop_count_tracker << "] "
                      << " Intra " << *running_intra_bucket_comm
                      << " full " << *running_insert_in_full
                      << " Total " << *running_intra_bucket_comm + *running_buffer_allocate + *running_local_compute + *running_all_to_all + *running_buffer_free + *running_insert_newt + *running_insert_in_full << std::endl;
        }
#endif

        batch_size--;
        loop_count_tracker++;

        if (iteration_count == 1)
            break;
    }

    delete[] offset;

    check_for_fixed_point(history);

    if (logging == true)
        print_all_relation();
}


#if 0
void RAM::execute_by_wall_clock(double batch_time, std::vector<u32>& history, std::map<u64, u64>& intern_map)
{
    double running_total = 0;
    double running_time = 0;

    while (running_time <= batch_time)
    {
        double elapsed_start = MPI_Wtime();

        double intra_start = MPI_Wtime();
        intra_bucket_comm_execute();
        double intra_end = MPI_Wtime();

        double allocate_buffers_start = MPI_Wtime();
        allocate_compute_buffers();
        double allocate_buffers_end = MPI_Wtime();

        double compute_start = MPI_Wtime();
        local_compute();
        double compute_end = MPI_Wtime();

        double all_to_all_start = MPI_Wtime();
        all_to_all();
        double all_to_all_end = MPI_Wtime();

        double free_buffers_start = MPI_Wtime();
        free_compute_buffers();
        double free_buffers_end = MPI_Wtime();

        double insert_in_newt_start = MPI_Wtime();
        local_insert_in_newt(intern_map);
        double insert_in_newt_end = MPI_Wtime();

        double insert_in_full_start = MPI_Wtime();
        local_insert_in_full();
        double insert_in_full_end = MPI_Wtime();


        running_total = running_total + (intra_end - intra_start) + (compute_end - compute_start) + (all_to_all_end - all_to_all_start) + (insert_in_newt_end - insert_in_newt_start) + (insert_in_full_end - insert_in_full_start);


        if (mcomm.get_rank() == 0)
            std::cout << "INNER [" << loop_count_tracker << "] "
                      << " Intra " << (intra_end - intra_start)
                      << " Buffer create " << (allocate_buffers_end - allocate_buffers_start)
                      << " compute " << (compute_end - compute_start)
                      << " All to all " << (all_to_all_end - all_to_all_start)
                      << " Buffer free " << (free_buffers_end - free_buffers_start)
                      << " Insert newt " << (insert_in_newt_end - insert_in_newt_start)
                      << " Insert full " << (insert_in_full_end - insert_in_full_start)
                      << " Total " << (intra_end - intra_start) + (compute_end - compute_start) + (all_to_all_end - all_to_all_start) + (insert_in_newt_end - insert_in_newt_start) + (insert_in_full_end - insert_in_full_start)
                      << " [ "
                      << running_time
                      << " ]" << std::endl;

        if (iteration_count == 1)
            break;

        double elapsed_end = MPI_Wtime();
        running_time = running_time + (elapsed_end - elapsed_start);

        double mint;
        MPI_Allreduce(&running_time, &mint, 1, MPI_DOUBLE, MPI_MAX, mcomm.get_local_comm());
        running_time = mint;

        loop_count_tracker++;
    }

    check_for_fixed_point(history);

    if (logging == true)
        print_all_relation();
}
#endif
