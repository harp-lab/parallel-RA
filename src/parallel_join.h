#ifndef PARALLEL_JOIN_H
#define PARALLEL_JOIN_H

#include <mpi.h>
#include <queue>
#include "google_btree_relation.h"
#include "balanced_hash_relation.h"
#include "balanced_parallel_io.h"





class parallel_RA
{

protected:
    int iteration_count = -1;
    u32 RA_type;
    mpi_comm mcomm;

public:

    parallel_RA ()
    {
    }

    virtual ~parallel_RA()
    {
    }

    void set_iteration_count(int ic)
    {
        iteration_count = ic;
        return;
    }

    int get_iteration_count()
    {
        return iteration_count;
    }

    void decrement_iteration_count()
    {
        iteration_count--;
    }

    virtual void set_join_input0(relation* i0, int g_type0) {return;}
    virtual relation* get_join_input0() {return NULL;}
    virtual int get_join_input0_graph_type() {return 0;}
    virtual void set_join_input1(relation* i1, int g_type1) {return;}
    virtual relation* get_join_input1() {return NULL;}
    virtual int get_join_input1_graph_type() {return 0;}
    virtual void set_join_output(relation*& out) {return;}
    virtual relation* get_join_output() {return NULL;}
    virtual void set_join_projection_index (int* projection_reorder_index_array, int projection_reorder_index_array_length) {return;}
    virtual void get_join_projection_index(int** projection_reorder_index_array, int* projection_reorder_index_array_length) {return;}
    virtual void set_join_column_count (int jcc) {return;}
    virtual int get_join_column_count () {return 0;}

    virtual void set_copy_input(relation* i0, int g_type0){return;}
    virtual relation* get_copy_input(){return NULL;}
    virtual int get_copy_input0_graph_type(){return 0;}
    virtual void set_copy_output(relation*& out) {return;}
    virtual relation* get_copy_output(){return NULL;}
    virtual void set_copy_rename_index (int* pria, int prial) {return;}
    virtual void get_copy_rename_index(int** projection_reorder_index_array, int* projection_reorder_index_array_length) {return;}


    void set_comm(mpi_comm& mcomm)
    {
        this->mcomm = mcomm;
    }


    u32 get_RA_type()
    {
        return RA_type;
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

    void intra_bucket_comm(u32 buckets,
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
            input_buffer[i] = vector_buffer_create_empty();

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
            vector_buffer_free(&input_buffer[i]);

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





#if 0
    u64 full_count(relation* output)
    {
        u64 sum = 0;
        //int full_element_count = output->get_full_element_count();
        u64 full_element_count = 0;
        int buckets = mcomm.get_number_of_buckets();
        u32* sub_bucket = output->get_sub_bucket_count();
        u32** sub_bucket_size = output->get_full_sub_bucket_element_count();

        for (int b = 0; b < buckets; b++)
        {
            for (u32 c = 0; c < sub_bucket[b]; c++)
                full_element_count = full_element_count + (u64)sub_bucket_size[b][c];
        }

        MPI_Allreduce(&full_element_count, &sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mcomm.get_comm());
        return sum;
    }
#endif


#if 0
    int print_full(relation* output)
    {
        int buckets = mcomm.get_number_of_buckets();
        int arity = output->get_arity();
        u32* bucket_map = output->get_bucket_map();
        google_relation *delta = output->get_delta();
        google_relation *full = output->get_full();
        u32 insert_success = 0;

        u64 t[arity];
        for (int i = 0; i < buckets; i++)
        {
            if (bucket_map[i] == 1)
            {
                for ( auto local_it = delta[i].begin(); local_it!= delta[i].end(); ++local_it )
                {
                    Map0* k = local_it->second;
                    for (auto it2 = k->begin(); it2 != k->end(); it2++)
                    {
                        t[0] = local_it->first;
                        t[1] = it2->first;

                        std::cout << t[0] << "\t" << t[1] << std::endl;
                    }
                }
            }
        }


        return insert_success;
    }
#endif



    void all_to_all(vector_buffer* local_join_output, int* process_size, u64 *outer_hash_buffer_size, u64 **outer_hash_data)
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
};

class parallel_join: public parallel_RA
{

private:

    relation* join_input0_table;
    int join_input0_graph_type;

    relation* join_input1_table;
    int join_input1_graph_type;

    relation* join_output_table;

    int join_column_count;

    std::vector<int> projection_reorder_index_array;
    int projection_reorder_index_array_length;

public:
    parallel_join()
    {
        RA_type = JOIN;
    }

    parallel_join(relation* G, int G_type, relation* T, int T_type, relation* output, int jc_count, std::vector<int> projection_reorder_index_array, int projection_reorder_index_array_length)
        : join_input0_table(G), join_input0_graph_type(G_type), join_input1_table(T), join_input1_graph_type(T_type), join_output_table(output), join_column_count(jc_count), projection_reorder_index_array(projection_reorder_index_array), projection_reorder_index_array_length(projection_reorder_index_array_length)
    {
        RA_type = JOIN;
    }


    void set_join_input0(relation* i0, int g_type0)
    {
        join_input0_table = i0;
        join_input0_graph_type = g_type0;
    }

    relation* get_join_input0()
    {
        return join_input0_table;
    }

    int get_join_input0_graph_type()
    {
        return join_input0_graph_type;
    }

    void set_join_input1(relation* i1, int g_type1)
    {
        join_input1_table = i1;
        join_input1_graph_type = g_type1;
    }

    relation* get_join_input1()
    {
        return join_input1_table;
    }

    int get_join_input1_graph_type()
    {
        return join_input1_graph_type;
    }

    void set_join_output(relation*& out)
    {
        join_output_table = out;
    }

    relation* get_join_output()
    {
        return join_output_table;
    }

    void set_join_projection_index (std::vector<int> pria, int prial)
    {
        projection_reorder_index_array_length = prial;
        projection_reorder_index_array = pria;
    }

    void set_join_column_count (int jcc)
    {
        join_column_count = jcc;
    }

    int get_join_column_count ()
    {
        return join_column_count;
    }

    void get_join_projection_index(std::vector<int>* projection_reorder_index_array, int* projection_reorder_index_array_length)
    {
        *projection_reorder_index_array_length = this->projection_reorder_index_array_length;
        *projection_reorder_index_array = this->projection_reorder_index_array;
    }

    // TODO join_order
    u32 local_join( u32 buckets,
                    int input0_buffer_size, int input0_buffer_width, u64 *input0_buffer, int join_order,
                    google_relation *input1, u32 i1_size, int input1_buffer_width,
                    int reorder_map_array_size, std::vector<int> reorder_map_array,
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

            vector_buffer temp_buffer = vector_buffer_create_empty();// vector_buffer_create_with_capacity(1024);
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
                    vector_buffer_append(&(local_join_output[iteration][index]), (const unsigned char*)projected_path, sizeof(u64)*input1_buffer_width);
                    local_join_inserts++;
                }
                else
                    local_join_duplicates++;
            }
            vector_buffer_free(&temp_buffer);
#endif

            sum1 = sum1 + (t2 - t1);
        }
        //std::cout << "local_join_inserts " << local_join_inserts << std::endl;
        deduplicate.remove_tuple();
        return local_join_duplicates;
    }
};

class parallel_copy: public parallel_RA
{

private:
    relation* copy_input0_table;
    int copy_input0_graph_type;
    relation* copy_output_table;

    std::vector<int> copy_reorder_index_array;
    int copy_reorder_index_array_length;

public:
    parallel_copy()
    {
        iteration_count = 1;
        RA_type = COPY;
    }

    parallel_copy(relation* G, int G_version, relation* output, std::vector<int> reorder_index_array, int reorder_index_array_size)
        : copy_input0_table(G), copy_input0_graph_type(G_version), copy_output_table(output), copy_reorder_index_array(reorder_index_array), copy_reorder_index_array_length(reorder_index_array_size)
    {
        RA_type = COPY;
    }



    void set_copy_input(relation* i0, int g_type0)
    {
        copy_input0_table = i0;
        copy_input0_graph_type = g_type0;
    }
    relation* get_copy_input()
    {
        return copy_input0_table;
    }
    int get_copy_input0_graph_type()
    {
        return copy_input0_graph_type;
    }
    void set_copy_output(relation*& out)
    {
        copy_output_table = out;
    }
    relation* get_copy_output()
    {
        return copy_output_table;
    }

    void set_copy_rename_index (std::vector<int> pria, int prial)
    {
        copy_reorder_index_array_length = prial;
        copy_reorder_index_array = pria;
    }

    void get_copy_rename_index(std::vector<int>* projection_reorder_index_array, int* projection_reorder_index_array_length)
    {
        *projection_reorder_index_array_length = this->copy_reorder_index_array_length;
        *projection_reorder_index_array = this->copy_reorder_index_array;
    }


    void local_copy(u32 buckets, google_relation* input, u32* input_bucket_map, relation* output, int reorder_map_length, std::vector<int> reorder_map, vector_buffer* local_join_output, int* process_size, int* cumulative_process_size)
    {
        u32 arity = output->get_arity();
        //assert(reorder_map_length == (int)arity);

        u32* output_sub_bucket_count = output->get_sub_bucket_per_bucket_count();
        u32** output_sub_bucket_rank = output->get_sub_bucket_rank();

        for (u32 i = 0; i < buckets; i++)
        {
            if (input_bucket_map[i] == 1)
            {
                vector_buffer temp_buffer = vector_buffer_create_empty();

                std::vector<u64> prefix = {};
                input[i].as_vector_buffer_recursive(&temp_buffer, prefix);

                //std::cout << "temp_buffer.size / sizeof(u64) " << temp_buffer.size / sizeof(u64) << std::endl;
                for (u32 s = 0; s < temp_buffer.size / sizeof(u64); s=s+arity)
                {
                    u64 reordered_cur_path[arity];
                    for (u32 j =0; j < arity; j++)
                        memcpy(reordered_cur_path + reorder_map[j], (&temp_buffer)->buffer + ((s + j) * sizeof(u64)), sizeof(u64));

                    uint64_t bucket_id = hash_function(reordered_cur_path[0]) % buckets;
                    uint64_t sub_bucket_id = hash_function(reordered_cur_path[1]) % output_sub_bucket_count[bucket_id];
                    int index = output_sub_bucket_rank[bucket_id][sub_bucket_id];

                    process_size[index] = process_size[index] + arity;
                    cumulative_process_size[index] = cumulative_process_size[index] + arity;

                    //if (mcomm.get_local_rank() == 1)
                    //    std::cout << s << " " << index << " " << reordered_cur_path[0] << " " << reordered_cur_path[1] << " " << process_size[index] << std::endl;

                    vector_buffer_append(&local_join_output[index], (const unsigned char*)reordered_cur_path, sizeof(u64)*arity);
                }
                vector_buffer_free(&temp_buffer);
            }
        }

        //std::cout << "Arity " << arity << std::endl;
        //std::cout << "PSSSSSSSSS " << process_size[0] << " " << process_size[1] << std::endl;

        return;
    }

};

#endif
