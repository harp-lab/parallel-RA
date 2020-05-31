#include "parallel_RA_inc.h"



u32 parallel_join::local_join( u32 buckets,
                int input0_buffer_size, int input0_buffer_width, u64 *input0_buffer, int join_order,
                google_relation *input1, u32 i1_size, int input1_buffer_width,
                std::vector<int> reorder_map_array,
                relation* output,
                vector_buffer** local_join_output, int** process_size, int* cumulative_process_size,
                u32 threshhold, int join_colun_count,
                u32* local_join_count, int iteration)
{
    u32 local_join_duplicates = 0;
    u32 local_join_inserts = 0;

    google_relation deduplicate;

    u32* output_sub_bucket_count = output->get_sub_bucket_per_bucket_count();
    u32** output_sub_bucket_rank = output->get_sub_bucket_rank();

    // TODO: fix this
    //int rank = mcomm.get_rank();


    //std::cout << "reorder_map_array_size " << reorder_map_array_size << std::endl;
    //std::cout << "input1_buffer_width " << input1_buffer_width << std::endl;
    //std::cout << "input0_buffer_width " << input0_buffer_width << std::endl;
    //std::cout << "join_colun_count " << join_colun_count << std::endl;
    //std::cout << "input0_buffer_size " << input0_buffer_size << std::endl;

    //assert(reorder_map_array_size == input1_buffer_width + input0_buffer_width - join_colun_count);
    //u64 reordered_cur_path[input1_buffer_width + input0_buffer_width - join_colun_count];
    u64 projected_path[input1_buffer_width];

    double t1, t2, sum1=0;
    //u64 sumc=0;


    for (int k1 = 0; k1 < input0_buffer_size; k1 = k1 + input0_buffer_width)
    {
        t1 = MPI_Wtime();
        u64 bucket_id = hash_function(input0_buffer[k1]) % buckets;


        std::vector<u64> prefix;
        prefix.reserve(1024);
        for (int jc=0; jc < join_colun_count; jc++)
            prefix.push_back(input0_buffer[k1 + jc]);

        vector_buffer temp_buffer;
        temp_buffer.vector_buffer_create_empty();// vector_buffer_create_with_capacity(1024);
        input1[bucket_id].as_vector_buffer_recursive(&temp_buffer, prefix);

        /*
        input1[bucket_id].as_vector_buffer_recursive_hack(prefix,
                                                          &local_join_inserts, &local_join_duplicates,
                                                          &deduplicate,
                                                          buckets, output_sub_bucket_count, output_sub_bucket_rank,
                                                          iteration,
                                                          local_join_output, process_size, cumulative_process_size,
                                                          input0_buffer[k1 + 1]);
        */

        //sumc = sumc + temp_buffer.size / sizeof(u64);
        t2 = MPI_Wtime();

#if 1
        for (u32 s = 0; s < temp_buffer.size / sizeof(u64); s = s + input1_buffer_width)
        {
            /*
            for (int i = 0; i < input1_buffer_width; i++)
                memcpy(reordered_cur_path + i, temp_buffer.buffer + (s + i)*sizeof(u64), sizeof(u64));

            for (int i = join_colun_count; i < input1_buffer_width; i++)
                reordered_cur_path[input1_buffer_width + (i - join_colun_count)] = input0_buffer[k1 + i];

            for (int i =0; i < input1_buffer_width + input0_buffer_width - join_colun_count; i++)
            {
                if (reorder_map_array[i] == -1)
                    continue;
                projected_path[reorder_map_array[i]] = reordered_cur_path[i];
            }
            */
            memcpy(projected_path, temp_buffer.buffer + (s+1)*sizeof(u64), sizeof(u64));
            projected_path[1] = input0_buffer[k1 + 1];

            //memcpy(projected_path + sizeof(u64), temp_buffer.buffer + (s+1)*sizeof(u64), sizeof(u64));
            //projected_path[0] = input0_buffer[k1 + 1];

            if (deduplicate.insert_tuple_from_array(projected_path, input1_buffer_width) == true)
            {
                uint64_t bucket_id = hash_function(projected_path[0]) % buckets;
                uint64_t sub_bucket_id = hash_function(projected_path[1]) % output_sub_bucket_count[bucket_id];
                int index = output_sub_bucket_rank[bucket_id][sub_bucket_id];

                process_size[iteration][index] = process_size[iteration][index] + input1_buffer_width;
                cumulative_process_size[index] = cumulative_process_size[index] + input1_buffer_width;
                local_join_output[iteration][index].vector_buffer_append((const unsigned char*)projected_path, sizeof(u64)*input1_buffer_width);
                local_join_inserts++;
            }
            else
                local_join_duplicates++;
        }
        temp_buffer.vector_buffer_free();
#endif

        sum1 = sum1 + (t2 - t1);
    }
    //std::cout << "local_join_inserts " << local_join_inserts << std::endl;
    deduplicate.remove_tuple();
    return local_join_duplicates;
}


/// Method implementing intra-bucket comm
///
/// Input:
///     rel (data that needs to be transmitted)
///     input_distinct_sub_bucket_rank_count (the number of processes a process has to send data to (for every bucket))
///     input_distinct_sub_bucket_rank (rank of the processes a process has to send data to (for every bucket))
///     input_bucket_map
///     output_distinct_sub_bucket_rank_count (the number of processes a process has to send data to (for every bucket))
///     output_distinct_sub_bucket_rank (rank of the processes a process has to send data to (for every bucket))
///     output_bucket_map
/// Output:
///     total_buffer_size
///     recvbuf

void parallel_join::intra_bucket_comm(u32 buckets,
                       google_relation *rel,
                       int* input_distinct_sub_bucket_rank_count, int** input_distinct_sub_bucket_rank, u32* input_bucket_map,
                       int* output_distinct_sub_bucket_rank_count, int** output_distinct_sub_bucket_rank, u32* output_bucket_map,
                       u64 *total_buffer_size, u64 **recvbuf)
{
    // buffer to hold relation data to be sent out
    vector_buffer *input_buffer = new vector_buffer[buckets];
    int *input_buffer_size = new int[buckets];

    //std::cout << "Buckets " << buckets << std::endl;

    u32** meta_buffer_size = new u32*[buckets];
    memset(meta_buffer_size, 0, sizeof(u32*) * buckets);

    *total_buffer_size = 0;
    u32* bucket_offset = new u32[buckets];


    u64 total_send_buffer_size = 0;
    for (u32 i = 0; i < buckets; i++)
    {
        // Buffer to store relation data
        input_buffer[i].vector_buffer_create_empty();

        // Puts btree data into a vector
        std::vector<u64> prefix = {};
        rel[i].as_vector_buffer_recursive(&(input_buffer[i]), prefix);

        // size of data to be sent
        input_buffer_size[i] = (&input_buffer[i])->size / sizeof(u64);
        total_send_buffer_size = total_send_buffer_size + input_buffer_size[i];


        meta_buffer_size[i] = new u32[input_distinct_sub_bucket_rank_count[i]];
        memset(meta_buffer_size[i], 0, sizeof(u32) * input_distinct_sub_bucket_rank_count[i]);

        u32 req_counter1 = 0;
        MPI_Request *req1 = new MPI_Request[output_distinct_sub_bucket_rank_count[i] + input_distinct_sub_bucket_rank_count[i]];
        MPI_Status *stat1 = new MPI_Status[output_distinct_sub_bucket_rank_count[i] + input_distinct_sub_bucket_rank_count[i]];

        if (input_bucket_map[i] == 1)
        {
            for (int r = 0; r < output_distinct_sub_bucket_rank_count[i]; r++)
            {
                int buffer_size = input_buffer_size[i];
                MPI_Isend(&buffer_size, 1, MPI_INT, output_distinct_sub_bucket_rank[i][r], 123, mcomm.get_local_comm(), &req1[req_counter1]);
                req_counter1++;
            }
        }

        if (output_bucket_map[i] == 1)
        {
            for (int r = 0; r < input_distinct_sub_bucket_rank_count[i]; r++)
            {
                MPI_Irecv(meta_buffer_size[i] + r, 1, MPI_INT, input_distinct_sub_bucket_rank[i][r], 123, mcomm.get_local_comm(), &req1[req_counter1]);
                req_counter1++;
            }
        }

        MPI_Waitall(req_counter1, req1, stat1);

        bucket_offset[i] = *total_buffer_size;
        for (int r = 0; r < input_distinct_sub_bucket_rank_count[i]; r++)
            *total_buffer_size = *total_buffer_size + meta_buffer_size[i][r];

        delete[] req1;
        delete[] stat1;
    }


    /*
    int rank = mcomm.get_rank();
    // Code to verify that the intra-bucket comm is setup correctly
    u64 global_send_buffer_size1 = 0;
    u64 global_send_buffer_size2 = 0;
    MPI_Allreduce(&total_send_buffer_size, &global_send_buffer_size1, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mcomm.get_local_comm());
    MPI_Allreduce(total_buffer_size, &global_send_buffer_size2, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mcomm.get_local_comm());
    //std::cout << "total_send_buffer_size = " << total_send_buffer_size << std::endl;
    //std::cout << "*total_buffer_size = " << *total_buffer_size << std::endl;

    assert(global_send_buffer_size1 == global_send_buffer_size2);
    if (rank == 0)
        std::cout << "[VERIFY intra_bucket] " << global_send_buffer_size1 << " " << global_send_buffer_size2 << std::endl;
    */


    /// Actual data Exchange

    // Allocate buffer
    *recvbuf = new u64[*total_buffer_size];

    // Non-blocking point-to-point data exchange
    for (u32 i = 0; i < buckets; i++)
    {
        u32 req_counter2 = 0;
        MPI_Request *req2 = new MPI_Request[output_distinct_sub_bucket_rank_count[i] + input_distinct_sub_bucket_rank_count[i]];
        MPI_Status *stat2 = new MPI_Status[output_distinct_sub_bucket_rank_count[i] + input_distinct_sub_bucket_rank_count[i]];

        // data
        if (input_bucket_map[i] == 1)
        {
            for (int r = 0; r < output_distinct_sub_bucket_rank_count[i]; r++)
            {
                if (input_buffer_size[i] != 0)
                {
                    MPI_Isend(input_buffer[i].buffer, input_buffer_size[i], MPI_UNSIGNED_LONG_LONG, output_distinct_sub_bucket_rank[i][r], 123, mcomm.get_local_comm(), &req2[req_counter2]);
                    req_counter2++;
                }
            }
        }

        u32 offset = 0;
        if (output_bucket_map[i] == 1)
        {
            for (int r = 0; r < input_distinct_sub_bucket_rank_count[i]; r++)
            {
                if (meta_buffer_size[i][r] != 0)
                {
                    MPI_Irecv((*recvbuf) + offset + bucket_offset[i], meta_buffer_size[i][r], MPI_UNSIGNED_LONG_LONG, input_distinct_sub_bucket_rank[i][r], 123, mcomm.get_local_comm(), &req2[req_counter2]);
                    offset = offset + meta_buffer_size[i][r];
                    req_counter2++;
                }
            }
        }

        MPI_Waitall(req_counter2, req2, stat2);
        input_buffer[i].vector_buffer_free();

        delete[] req2;
        delete[] stat2;
        delete[] meta_buffer_size[i];
    }
    delete[] meta_buffer_size;
    delete[] input_buffer;
    delete[] input_buffer_size;
    delete[] bucket_offset;

    return;

}
