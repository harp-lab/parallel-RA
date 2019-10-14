#ifndef PARALLEL_JOIN_H
#define PARALLEL_JOIN_H

#include <mpi.h>
#include <queue>
#include "balanced_hash_relation.h"
#include "balanced_parallel_io.h"
#include "btree/btree_map.h"



enum {STATIC, DYNAMIC};

class parallel_RA
{
private:
    u32 RA_type;
    //int projection_index[3];
    int projection_index_a;
    int projection_index_b;
    int projection_index_c;

    /* JOIN */
    relation* join_input0_table;
    u32 join_input0_graph_type;

    relation* join_input1_table;
    u32 join_input1_graph_type;

    relation* join_output_table;


    /* COPY */
    relation* copy_input0_table;
    u32 copy_input0_graph_type;

    relation* copy_output_table;


    mpi_comm mcomm;


public:

    parallel_RA (int RA_Type, mpi_comm mcomm)
    {
        this->RA_type = RA_Type;
        this->mcomm = mcomm;
    }

    void get_join_projection_index(int* pi1, int* pi2, int* pi3)
    {
        *pi1 = projection_index_a;
        *pi2 = projection_index_b;
        *pi3 = projection_index_c;

        return;
    }

    void get_copy_projection_index(int* pi1, int* pi2)
    {
        *pi1 = projection_index_a;
        *pi2 = projection_index_b;

        return;
    }

    // Join
    void join_input0(relation* i0, int g_type0)
    {
        //join_input0_table = new relation();
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

    void join_input1(relation* i1, int g_type1)
    {
        //join_input1_table = new relation();
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

    void join_output(relation* out)
    {
        //join_output_table = new relation();
        join_output_table = out;
    }
    relation* get_join_output()
    {
        return join_output_table;
    }


    // Copy
    void copy_input(relation* i0, int g_type0)
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
    void copy_output(relation* out)
    {
        copy_output_table = out;
    }
    relation* get_copy_output()
    {
        return copy_output_table;
    }


    void set_join_projection_index (int pi1, int pi2, int pi3)
    {
        this->projection_index_a = pi1;
        this->projection_index_b = pi2;
        this->projection_index_c = pi3;
    }

    void set_copy_projection_index (int pi1, int pi2)
    {
        this->projection_index_a = pi1;
        this->projection_index_b = pi2;
    }

    u32 get_RA_type() {return RA_type;}



    bool fixed_point_check(relation* output)
    {
        bool fixed_point1 = false;
        int sum = 0;
        int delta_element_count = output->get_delta_element_count();
        MPI_Allreduce(&delta_element_count, &sum, 1, MPI_INT, MPI_SUM, mcomm.get_comm());

        bool fixed_point2 = false;
        int sum1 = 0;
        int full_inserts = output->get_full_inserts_element_count();
        MPI_Allreduce(&full_inserts, &sum1, 1, MPI_INT, MPI_SUM, mcomm.get_comm());

        if(sum == 0)
            fixed_point1 = true;

        if(sum1 == 0)
            fixed_point2 = true;

        return fixed_point1 & fixed_point2;
    }


    u32 full_count(relation* output)
    {
        int sum = 0;
        //int full_element_count = output->get_full_element_count();
        int full_element_count = 0;
        int buckets = mcomm.get_number_of_buckets();
        u32* sub_bucket = output->get_sub_bucket_count();
        u32** sub_bucket_size = output->get_full_sub_bucket_element_count();

        for (int b = 0; b < buckets; b++)
        {
            for (u32 c = 0; c < sub_bucket[b]; c++)
                full_element_count = full_element_count + sub_bucket_size[b][c];
        }

        MPI_Allreduce(&full_element_count, &sum, 1, MPI_INT, MPI_SUM, mcomm.get_comm());
        return sum;
    }



    void clique_comm(int input_arity, google_relation *full, u32 full_size, int* input_distinct_sub_bucket_rank_count, int** input_distinct_sub_bucket_rank, u32* input_bucket_map, relation* output, u64 *total_buffer_size, u64 **recvbuf)
    {
#if 1
        u32 buckets = mcomm.get_number_of_buckets();
        //int rank = mcomm.get_rank();
        //u32 nprocs = mcomm.get_nprocs();

        vector_buffer *input_buffer = new vector_buffer[buckets];
        int *input_buffer_size = new int[buckets];

        int* output_distinct_sub_bucket_rank_count;
        int** output_distinct_sub_bucket_rank;
        u32* output_bucket_map;

        u32** meta_buffer_size = new u32*[buckets];
        memset(meta_buffer_size, 0, sizeof(u32*) * buckets);

        *total_buffer_size = 0;
        u32* bucket_offset = new u32[buckets];

        u64 val[input_arity];

        output_distinct_sub_bucket_rank_count = output->get_distinct_sub_bucket_rank_count();
        output_distinct_sub_bucket_rank = output->get_distinct_sub_bucket_rank();
        output_bucket_map = output->get_bucket_map();

        u32 total_send_buffer_size = 0;
        for (u32 i = 0; i < buckets; i++)
        {
            input_buffer[i] = vector_buffer_create_empty();
            for (auto it = full[i].begin(); it != full[i].end(); it++)
            {
                Map0* it2 = it->second;
                for (auto dit2 = it2->begin(); dit2 != it2->end(); dit2++)
                {
                    val[0] = it->first;
                    val[1] = dit2->first;
                    vector_buffer_append(&input_buffer[i],(unsigned char *) val, sizeof(u64)*input_arity);
                }
            }

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
                    MPI_Isend(&buffer_size, 1, MPI_INT, output_distinct_sub_bucket_rank[i][r], 123, MPI_COMM_WORLD, &req1[req_counter1]);
                    req_counter1++;
                }
            }

            if (output_bucket_map[i] == 1)
            {
                for (int r = 0; r < input_distinct_sub_bucket_rank_count[i]; r++)
                {
                    MPI_Irecv(meta_buffer_size[i] + r, 1, MPI_INT, input_distinct_sub_bucket_rank[i][r], 123, MPI_COMM_WORLD, &req1[req_counter1]);
                    req_counter1++;
                }
            }

            MPI_Waitall(req_counter1, req1, stat1);

            bucket_offset[i] = *total_buffer_size;
            for (int r = 0; r < input_distinct_sub_bucket_rank_count[i]; r++)
            {
                //if (rank == 0)
                //    std::cout << "meta_buffer_size[i][r] " << meta_buffer_size[i][r] << std::endl;
                *total_buffer_size = *total_buffer_size + meta_buffer_size[i][r];
            }

            delete[] req1;
            delete[] stat1;
        }

#if 0
        u64 global_send_buffer_size1 = 0;
        u64 global_send_buffer_size2 = 0;
        MPI_Allreduce(&total_send_buffer_size, &global_send_buffer_size1, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(total_buffer_size, &global_send_buffer_size2, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        std::cout << "total_send_buffer_size = " << total_send_buffer_size << std::endl;
        std::cout << "*total_buffer_size = " << *total_buffer_size << std::endl;
        if (rank == 0)
        {
            std::cout << "[VERIFY CLIQUE] " << global_send_buffer_size1 << " " << global_send_buffer_size2 << std::endl;
        }
#endif

        /* Allocate buffer */
        //std::cout << "CLIQUE BUFFER SIZE " << *total_buffer_size << std::endl;
        *recvbuf = new u64[*total_buffer_size];

        /* Actual data exchange */
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
                        MPI_Isend(input_buffer[i].buffer, input_buffer_size[i], MPI_UNSIGNED_LONG_LONG, output_distinct_sub_bucket_rank[i][r], 123, MPI_COMM_WORLD, &req2[req_counter2]);
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
                        MPI_Irecv((*recvbuf) + offset + bucket_offset[i], meta_buffer_size[i][r], MPI_UNSIGNED_LONG_LONG, input_distinct_sub_bucket_rank[i][r], 123, MPI_COMM_WORLD, &req2[req_counter2]);
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
#endif
#if 0
        *recvbuf = new u64[full_size * 2];
        *total_buffer_size = full_size * 2;
        u32 c = 0;
        u64 val[2];
        //for (u32 i = 0; i < buckets; i++)
        u32 i = mcomm.get_rank();
        {
            for (auto it = full[i].begin(); it != full[i].end(); it++)
            {
                Map0* it2 = it->second;
                for (auto dit2 = it2->begin(); dit2 != it2->end(); dit2++)
                {
                    val[0] = it->first;
                    val[1] = dit2->first;
                    (*recvbuf)[c] = val[0];
                    (*recvbuf)[c+1] = val[1];
                    c = c+2;
                }
            }
        }
        //std::cout << "full_size " << full_size * 2<< std::endl;
        //std::cout << "c " << c<< std::endl;
#endif
        return;
    }

    void local_copy(relation* input, relation* output, vector_buffer* local_join_output, int* process_size)
    {
        int rank = mcomm.get_rank();
        u32 output_arity = output->get_arity();
        u32* output_sub_bucket_count = output->get_sub_bucket_count();
        u32** output_sub_bucket_rank = output->get_sub_bucket_rank();
        int buckets = mcomm.get_number_of_buckets();
        google_relation *delta = input->get_delta();
        u32* bucket_map = input->get_bucket_map();
        u64 val[2] = {0, 0};

        u32 input0_buffer_size = 0;


        //int i = rank;
        //for (int k1 = 0; k1 < input0_buffer_size; k1=k1+input0_arity)
        for (int i = 0; i < buckets; i++)
        {
            if (bucket_map[i] == 1)
            {
                for ( auto local_it = delta[i].begin(); local_it!= delta[i].end(); ++local_it )
                {
                    Map0* k = local_it->second;
                    for (auto it2 = k->begin(); it2 != k->end(); it2++)
                    {
                        val[0] = it2->first;
                        val[1] = local_it->first;

                        uint64_t bucket_id = hash_function(val[0]) % buckets;
                        uint64_t sub_bucket_id = hash_function(val[1]) % output_sub_bucket_count[bucket_id];

                        int index = output_sub_bucket_rank[bucket_id][sub_bucket_id];
                        vector_buffer_append(&local_join_output[index], (unsigned char *) val, sizeof(u64) * output_arity);
                        process_size[index] = process_size[index] + output_arity;
                        input0_buffer_size++;
                    }
                }
            }
        }

        if (rank == 0)
        std::cout << "[Local Copy] " << input0_buffer_size << std::endl;

        return;
    }

    void cumulative_local_copy(relation* input, relation* output, vector_buffer** local_join_output, int** process_size, int* cumulative_process_size, int iteration)
    {
        int rank = mcomm.get_rank();
        u32 output_arity = output->get_arity();
        u32* output_sub_bucket_count = output->get_sub_bucket_count();
        u32** output_sub_bucket_rank = output->get_sub_bucket_rank();
        int buckets = mcomm.get_number_of_buckets();
        google_relation *delta = input->get_delta();
        u32* bucket_map = input->get_bucket_map();
        u64 val[2] = {0, 0};

        u32 input0_buffer_size = 0;


        //int cy1 = 0;
        //int i = rank;
        //for (int k1 = 0; k1 < input0_buffer_size; k1=k1+input0_arity)
        for (int i = 0; i < buckets; i++)
        {
            if (bucket_map[i] == 1)
            {
                for ( auto local_it = delta[i].begin(); local_it!= delta[i].end(); ++local_it )
                {
                    Map0* k = local_it->second;
                    for (auto it2 = k->begin(); it2 != k->end(); it2++)
                    {
                        val[0] = it2->first;
                        val[1] = local_it->first;

                        uint64_t bucket_id = hash_function(val[0]) % buckets;
                        uint64_t sub_bucket_id = hash_function(val[1]) % output_sub_bucket_count[bucket_id];

                        int index = output_sub_bucket_rank[bucket_id][sub_bucket_id];
                        vector_buffer_append(&local_join_output[iteration][index], (unsigned char *) val, sizeof(u64) * output_arity);
                        process_size[iteration][index] = process_size[iteration][index] + output_arity;
                        cumulative_process_size[index] = cumulative_process_size[index] + output_arity;
                        //std::cout << "C " << cy1 << std::endl;
                        //cy1++;
                        input0_buffer_size++;
                    }
                }
            }
        }

        if (rank == 0)
        std::cout << "[Local Copy] " << input0_buffer_size << std::endl;

        return;
    }


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



    bool local_join(u32 threshhold, int input0_buffer_size, int *offset, int input0_arity, u64 *input0_buffer, google_relation *input1, u32 i1_size, relation* output, vector_buffer* local_join_output, int* process_size, int projection_and_rename_b, int projection_and_rename_c, u32* local_join_count, int iteration)
    {
        u32 local_join_duplicates = 0;
        u32 local_join_inserts = 0;
        google_relation tempT;
        u32 output_arity = output->get_arity();
        u32* output_sub_bucket_count = output->get_sub_bucket_count();
        u32** output_sub_bucket_rank = output->get_sub_bucket_rank();
        int buckets = mcomm.get_number_of_buckets();
        u64 val[2] = {0, 0};
        int elements_accessed = 0;
        int rank = mcomm.get_rank();

        //std::cout << "Entering JOIN " << rank << std::endl;

        if (*offset > input0_buffer_size || input0_buffer_size == 0 || i1_size == 0)
        {
            if (rank == 0)
                std::cout  <<"[Join Done] [Local Join] " << i1_size << " " << elements_accessed << " (" << input0_buffer_size << ") " << local_join_inserts << " Duplicates " << local_join_duplicates << " Offset " << *offset << std::endl;

            return true;
        }


        /*
        if (iteration == 1 && rank == 0)
            for (int k1 = *offset; k1 < input0_buffer_size; k1 = k1 + input0_arity)
                std::cout <<"T " << input0_buffer[k1] << " " << input0_buffer[k1 + 1] << std::endl;

        u64 t[2];
        for (int i = 0; i < buckets; i++)
        {
            //if (bucket_map[i] == 1)
            {
                for ( auto local_it = input1[i].begin(); local_it!= input1[i].end(); ++local_it )
                {
                    Map0* k = local_it->second;
                    for (auto it2 = k->begin(); it2 != k->end(); it2++)
                    {
                        t[0] = local_it->first;
                        t[1] = it2->first;

                        if (iteration == 1 && rank == 0)
                            std::cout << "G " << t[0] << " " << t[1] << std::endl;
                    }
                }
            }
        }
        */



        if (projection_and_rename_b == 1 && projection_and_rename_c == 0)
        {
            for (int k1 = *offset; k1 < input0_buffer_size; k1 = k1 + input0_arity)
            {
                //if (iteration == 1 && rank == 0)
                    //std::cout << "Input " << input0_buffer[k1] << " " << input0_buffer[k1 + 1] << std::endl;
                int i = hash_function(input0_buffer[k1]) % buckets;
                //for (int i = 0; i < buckets; i++)
                {
                    auto itd = input1[i].find(input0_buffer[k1]);
                    if( itd != input1[i].end() )
                    {
                        Map0* Git = itd->second;
                        for (auto it2 = Git->begin(); it2 != Git->end(); it2++)
                        {
                            auto itx = tempT.find(input0_buffer[k1 + 1]);
                            if( itx != tempT.end() )
                            {
                                auto it2x = (itx->second)->find(it2->first);
                                if( it2x != (itx->second)->end() )
                                {
                                    local_join_duplicates++;
                                }
                                else
                                {
                                    (itx->second)->insert(std::make_pair(it2->first, 0));
                                    tempT[input0_buffer[k1 + 1]] = itx->second;
                                    local_join_inserts++;

                                    val[1] = it2->first;
                                    val[0] = input0_buffer[k1 + 1];
                                    //if (iteration == 1)
                                        //if (iteration == 1 && rank == 0)
                                        //std::cout << rank << " Matches " << input0_buffer[k1] << " " << input0_buffer[k1 + 1] << " | " << input0_buffer[k1] << " " <<  it2->first << std::endl;
                                    uint64_t bucket_id = hash_function(val[0]) % buckets;
                                    uint64_t sub_bucket_id = hash_function(val[1]) % output_sub_bucket_count[bucket_id];
                                    int index = output_sub_bucket_rank[bucket_id][sub_bucket_id];
                                    vector_buffer_append(&(local_join_output[index]), (unsigned char *) val, sizeof(u64) * output_arity);
                                    process_size[index] = process_size[index] + output_arity;
                                }
                            }
                            else
                            {
                                Map0* k = new Map0();
                                k->insert(std::make_pair(it2->first, 0));
                                tempT[input0_buffer[k1 + 1]] = k;
                                local_join_inserts++;

                                val[1] = it2->first;
                                val[0] = input0_buffer[k1 + 1];
                                //if (iteration == 1)
                                    //if (iteration == 1 && rank == 0)
                                    //std::cout << rank << " Matches " << input0_buffer[k1] << " " << input0_buffer[k1 + 1] << " | " << input0_buffer[k1] << " " <<  it2->first << std::endl;
                                uint64_t bucket_id = hash_function(val[0]) % buckets;
                                uint64_t sub_bucket_id = hash_function(val[1]) % output_sub_bucket_count[bucket_id];
                                int index = output_sub_bucket_rank[bucket_id][sub_bucket_id];
                                vector_buffer_append(&(local_join_output[index]), (unsigned char *) val, sizeof(u64) * output_arity);
                                process_size[index] = process_size[index] + output_arity;

                            }
                        }
                    }
                }

                if (local_join_inserts > threshhold)
                {
                    //std::cout << "FINAL1 " << input0_buffer[k1] << " " << input0_buffer[k1 + 1] << " " << k1 << " Temp Count " << temp_count << std::endl;

                    if (rank == 0)
                        std::cout << "[Threshold reached A] [Local Join] " << i1_size << " " << elements_accessed << " (" << input0_buffer_size << ") " << local_join_inserts << " Offset "<< *offset << " Duplicates " << local_join_duplicates << std::endl;

                    for(google_relation::iterator ix = tempT.begin(); ix != tempT.end(); ix++)
                        delete (ix->second);

                    *offset = k1 + input0_arity;
                    *local_join_count = local_join_inserts;
                    return false;
                }
                elements_accessed++;

                //std::cout << "FINAL2 " << input0_buffer[k1] << " " << input0_buffer[k1 + 1] << " " << k1 << " Temp Count " << temp_count << std::endl;

            }
        }
        else if (projection_and_rename_b == 0 && projection_and_rename_c == 1)
        {
            for (int k1 = *offset; k1 < input0_buffer_size; k1 = k1 + input0_arity)
            {
                //if (iteration == 1 && rank == 0)
                    //std::cout << "Input " << input0_buffer[k1] << " " << input0_buffer[k1 + 1] << std::endl;
                //int i = rank;
                //for (int i = 0; i < buckets; i++)
                int i = hash_function(input0_buffer[k1]) % buckets;
                {
                    auto itd = input1[i].find(input0_buffer[k1]);
                    if( itd != input1[i].end() )
                    {
                        Map0* Git = itd->second;
                        for (auto it2 = Git->begin(); it2 != Git->end(); it2++)
                        {
                            auto itx = tempT.find(it2->first);
                            if( itx != tempT.end() )
                            {
                                auto it2x = (itx->second)->find(input0_buffer[k1 + 1]);
                                if( it2x != (itx->second)->end() )
                                {
                                    local_join_duplicates++;
                                }
                                else
                                {
                                    (itx->second)->insert(std::make_pair(input0_buffer[k1 + 1], 0));
                                    tempT[it2->first] = itx->second;
                                    local_join_inserts++;

                                    val[0] = it2->first;
                                    val[1] = input0_buffer[k1 + 1];
                                    //if (iteration == 1)
                                        //if (iteration == 1 && rank == 0)
                                        //std::cout << rank << " Matches " << input0_buffer[k1] << " " << input0_buffer[k1 + 1] << " | " << input0_buffer[k1] << " " <<  it2->first << std::endl;
                                    uint64_t bucket_id = hash_function(val[0]) % buckets;
                                    uint64_t sub_bucket_id = hash_function(val[1]) % output_sub_bucket_count[bucket_id];
                                    int index = output_sub_bucket_rank[bucket_id][sub_bucket_id];
                                    vector_buffer_append(&(local_join_output[index]), (unsigned char *) val, sizeof(u64) * output_arity);
                                    process_size[index] = process_size[index] + output_arity;

                                }
                            }
                            else
                            {
                                Map0* k = new Map0();
                                k->insert(std::make_pair(input0_buffer[k1 + 1], 0));
                                tempT[it2->first] = k;
                                local_join_inserts++;

                                val[0] = it2->first;
                                val[1] = input0_buffer[k1 + 1];
                                //if (iteration == 1)
                                    //if (iteration == 1 && rank == 0)
                                    //std::cout << rank << " Matches " << input0_buffer[k1] << " " << input0_buffer[k1 + 1] << " | " << input0_buffer[k1] << " " <<  it2->first << std::endl;
                                uint64_t bucket_id = hash_function(val[0]) % buckets;
                                uint64_t sub_bucket_id = hash_function(val[1]) % output_sub_bucket_count[bucket_id];
                                int index = output_sub_bucket_rank[bucket_id][sub_bucket_id];
                                vector_buffer_append(&(local_join_output[index]), (unsigned char *) val, sizeof(u64) * output_arity);
                                process_size[index] = process_size[index] + output_arity;

                            }
                        }
                    }
                }

                if (local_join_inserts > threshhold)
                {
                    //std::cout << "FINAL1 " << input0_buffer[k1] << " " << input0_buffer[k1 + 1] << " " << k1 << " Temp Count " << temp_count << std::endl;

                    if (rank == 0)
                        std::cout << "[Threshold reached B] [Local Join] " << i1_size << " " << elements_accessed << " (" << input0_buffer_size << ") " << local_join_inserts << " Offset "<< *offset << " Duplicates " << local_join_duplicates << std::endl;

                    for(google_relation::iterator ix = tempT.begin(); ix != tempT.end(); ix++)
                        delete (ix->second);

                    *offset = k1 + input0_arity;
                    *local_join_count = local_join_inserts;
                    return false;
                }
                elements_accessed++;

                //std::cout << "FINAL2 " << input0_buffer[k1] << " " << input0_buffer[k1 + 1] << " " << k1 << " Temp Count " << temp_count << std::endl;

            }

        }
        if (rank == 0)
            std::cout  <<"[Join Complete] [Local Join] " << i1_size << " " << elements_accessed << " (" << input0_buffer_size << ") " << local_join_inserts << " Duplicates " << local_join_duplicates << " Offset " << *offset << std::endl;

        for(google_relation::iterator ix = tempT.begin(); ix != tempT.end(); ix++)
            delete (ix->second);

        *local_join_count = local_join_inserts;
        //*offset = input0_buffer_size + 1;

        return true;
    }


    bool cumulative_local_join(u32 threshhold, int input0_buffer_size, int *offset, int input0_arity, u64 *input0_buffer, google_relation *input1, u32 i1_size, relation* output, vector_buffer** local_join_output, int** process_size, int* cumulative_process_size, int projection_and_rename_b, int projection_and_rename_c, u32* local_join_count, int iteration)
    {
        u32 local_join_duplicates = 0;
        u32 local_join_inserts = 0;
        google_relation tempT;
        u32 output_arity = output->get_arity();
        u32* output_sub_bucket_count = output->get_sub_bucket_count();
        u32** output_sub_bucket_rank = output->get_sub_bucket_rank();
        int buckets = mcomm.get_number_of_buckets();
        u64 val[2] = {0, 0};
        int elements_accessed = 0;
        int rank = mcomm.get_rank();

        if (*offset > input0_buffer_size || input0_buffer_size == 0 || i1_size == 0)
        {
            if (rank == 0)
                std::cout  <<"[Join Done] [Local Join] " << i1_size << " " << elements_accessed << " (" << input0_buffer_size << ") " << local_join_inserts << " Duplicates " << local_join_duplicates << " Offset " << *offset << std::endl;

            return true;
        }

        //int cx1 = 0;
        if (projection_and_rename_b == 1 && projection_and_rename_c == 0)
        {
            for (int k1 = *offset; k1 < input0_buffer_size; k1 = k1 + input0_arity)
            {
                //int i = rank;
                for (int i = 0; i < buckets; i++)
                {
                    auto itd = input1[i].find(input0_buffer[k1]);
                    if( itd != input1[i].end() )
                    {
                        Map0* Git = itd->second;
                        for (auto it2 = Git->begin(); it2 != Git->end(); it2++)
                        {
                            auto itx = tempT.find(input0_buffer[k1 + 1]);
                            if( itx != tempT.end() )
                            {
                                auto it2x = (itx->second)->find(it2->first);
                                if( it2x != (itx->second)->end() )
                                {
                                    local_join_duplicates++;
                                }
                                else
                                {
                                    (itx->second)->insert(std::make_pair(it2->first, 0));
                                    tempT[input0_buffer[k1 + 1]] = itx->second;
                                    local_join_inserts++;

                                    val[1] = it2->first;
                                    val[0] = input0_buffer[k1 + 1];
                                    //if (iteration == 1)
                                        //if (iteration == 1 && rank == 0)
                                        //std::cout << rank << " Matches " << input0_buffer[k1] << " " << input0_buffer[k1 + 1] << " | " << input0_buffer[k1] << " " <<  it2->first << std::endl;
                                    uint64_t bucket_id = hash_function(val[0]) % buckets;
                                    uint64_t sub_bucket_id = hash_function(val[1]) % output_sub_bucket_count[bucket_id];
                                    int index = output_sub_bucket_rank[bucket_id][sub_bucket_id];
                                    vector_buffer_append(&(local_join_output[iteration][index]), (unsigned char *) val, sizeof(u64) * output_arity);
                                    //std::cout << cx1 << std::endl;
                                    //cx1++;
                                    process_size[iteration][index] = process_size[iteration][index] + output_arity;
                                    cumulative_process_size[index] = cumulative_process_size[index] + output_arity;
                                }
                            }
                            else
                            {
                                Map0* k = new Map0();
                                k->insert(std::make_pair(it2->first, 0));
                                tempT[input0_buffer[k1 + 1]] = k;
                                local_join_inserts++;

                                val[1] = it2->first;
                                val[0] = input0_buffer[k1 + 1];
                                //if (iteration == 1)
                                    //if (iteration == 1 && rank == 0)
                                    //std::cout << rank << " Matches " << input0_buffer[k1] << " " << input0_buffer[k1 + 1] << " | " << input0_buffer[k1] << " " <<  it2->first << std::endl;
                                uint64_t bucket_id = hash_function(val[0]) % buckets;
                                uint64_t sub_bucket_id = hash_function(val[1]) % output_sub_bucket_count[bucket_id];
                                int index = output_sub_bucket_rank[bucket_id][sub_bucket_id];
                                vector_buffer_append(&(local_join_output[iteration][index]), (unsigned char *) val, sizeof(u64) * output_arity);
                                //std::cout << cx1 << std::endl;
                                //cx1++;
                                process_size[iteration][index] = process_size[iteration][index] + output_arity;
                                cumulative_process_size[index] = cumulative_process_size[index] + output_arity;
                            }
                        }
                    }
                }

                if (local_join_inserts > threshhold)
                {
                    //std::cout << "FINAL1 " << input0_buffer[k1] << " " << input0_buffer[k1 + 1] << " " << k1 << " Temp Count " << temp_count << std::endl;

                    if (rank == 0)
                        std::cout << "[Threshold reached A] [Local Join] " << i1_size << " " << elements_accessed << " (" << input0_buffer_size << ") " << local_join_inserts << " Offset "<< *offset << " Duplicates " << local_join_duplicates << std::endl;

                    for(google_relation::iterator ix = tempT.begin(); ix != tempT.end(); ix++)
                        delete (ix->second);

                    *offset = k1 + input0_arity;
                    *local_join_count = local_join_inserts;
                    return false;
                }
                elements_accessed++;

                //std::cout << "FINAL2 " << input0_buffer[k1] << " " << input0_buffer[k1 + 1] << " " << k1 << " Temp Count " << temp_count << std::endl;

            }
        }
        else if (projection_and_rename_b == 0 && projection_and_rename_c == 1)
        {
            for (int k1 = *offset; k1 < input0_buffer_size; k1 = k1 + input0_arity)
            {
                //if (iteration == 1 && rank == 0)
                    //std::cout << "Input " << input0_buffer[k1] << " " << input0_buffer[k1 + 1] << std::endl;
                //int i = rank;
                for (int i = 0; i < buckets; i++)
                {
                    auto itd = input1[i].find(input0_buffer[k1]);
                    if( itd != input1[i].end() )
                    {
                        Map0* Git = itd->second;
                        for (auto it2 = Git->begin(); it2 != Git->end(); it2++)
                        {
                            auto itx = tempT.find(it2->first);
                            if( itx != tempT.end() )
                            {
                                auto it2x = (itx->second)->find(input0_buffer[k1 + 1]);
                                if( it2x != (itx->second)->end() )
                                {
                                    local_join_duplicates++;
                                }
                                else
                                {
                                    (itx->second)->insert(std::make_pair(input0_buffer[k1 + 1], 0));
                                    tempT[it2->first] = itx->second;
                                    local_join_inserts++;

                                    val[0] = it2->first;
                                    val[1] = input0_buffer[k1 + 1];
                                    //if (iteration == 1)
                                        //if (iteration == 1 && rank == 0)
                                        //std::cout << rank << " Matches " << input0_buffer[k1] << " " << input0_buffer[k1 + 1] << " | " << input0_buffer[k1] << " " <<  it2->first << std::endl;
                                    uint64_t bucket_id = hash_function(val[0]) % buckets;
                                    uint64_t sub_bucket_id = hash_function(val[1]) % output_sub_bucket_count[bucket_id];
                                    int index = output_sub_bucket_rank[bucket_id][sub_bucket_id];
                                    vector_buffer_append(&(local_join_output[iteration][index]), (unsigned char *) val, sizeof(u64) * output_arity);
                                    //std::cout << cx1 << std::endl;
                                    //cx1++;
                                    process_size[iteration][index] = process_size[iteration][index] + output_arity;
                                    cumulative_process_size[index] = cumulative_process_size[index] + output_arity;
                                }
                            }
                            else
                            {
                                Map0* k = new Map0();
                                k->insert(std::make_pair(input0_buffer[k1 + 1], 0));
                                tempT[it2->first] = k;
                                local_join_inserts++;

                                val[0] = it2->first;
                                val[1] = input0_buffer[k1 + 1];
                                //if (iteration == 1)
                                    //if (iteration == 1 && rank == 0)
                                    //std::cout << rank << " Matches " << input0_buffer[k1] << " " << input0_buffer[k1 + 1] << " | " << input0_buffer[k1] << " " <<  it2->first << std::endl;
                                uint64_t bucket_id = hash_function(val[0]) % buckets;
                                uint64_t sub_bucket_id = hash_function(val[1]) % output_sub_bucket_count[bucket_id];
                                int index = output_sub_bucket_rank[bucket_id][sub_bucket_id];
                                vector_buffer_append(&(local_join_output[iteration][index]), (unsigned char *) val, sizeof(u64) * output_arity);
                                //std::cout << cx1 << std::endl;
                                //cx1++;
                                process_size[iteration][index] = process_size[iteration][index] + output_arity;
                                cumulative_process_size[index] = cumulative_process_size[index] + output_arity;
                            }
                        }
                    }
                }

                if (local_join_inserts > threshhold)
                {
                    //std::cout << "FINAL1 " << input0_buffer[k1] << " " << input0_buffer[k1 + 1] << " " << k1 << " Temp Count " << temp_count << std::endl;

                    if (rank == 0)
                        std::cout << "[Threshold reached B] [Local Join] " << i1_size << " " << elements_accessed << " (" << input0_buffer_size << ") " << local_join_inserts << " Offset "<< *offset << " Duplicates " << local_join_duplicates << std::endl;

                    for(google_relation::iterator ix = tempT.begin(); ix != tempT.end(); ix++)
                        delete (ix->second);

                    *offset = k1 + input0_arity;
                    *local_join_count = local_join_inserts;
                    return false;
                }
                elements_accessed++;

                //std::cout << "FINAL2 " << input0_buffer[k1] << " " << input0_buffer[k1 + 1] << " " << k1 << " Temp Count " << temp_count << std::endl;

            }

        }
        if (rank == 0)
            std::cout  <<"[Join Complete] [Local Join] " << i1_size << " " << elements_accessed << " (" << input0_buffer_size << ") " << local_join_inserts << " Duplicates " << local_join_duplicates << " Offset " << *offset << " Duplicates " << local_join_duplicates << std::endl;

        for(google_relation::iterator ix = tempT.begin(); ix != tempT.end(); ix++)
            delete (ix->second);

        *local_join_count = local_join_inserts;
        *offset = input0_buffer_size + 1;

        return true;
    }





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




    int insert_in_newt(u64 outer_hash_buffer_size, u64 *outer_hash_data, relation* output)
    {
        int rank = mcomm.get_rank();
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

        if (rank == 0)
        {
            int new_count = output->get_new_element_count();
            std::cout << "[A] Local Inserts in new " << successful_insert << " Tried to inset " << outer_hash_buffer_size/2 << " New count (" << new_count << ")" << std::endl;
        }


        return successful_insert;
    }






};
#endif
