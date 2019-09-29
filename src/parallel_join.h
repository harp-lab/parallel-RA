#ifndef PARALLEL_JOIN_H
#define PARALLEL_JOIN_H

#include <mpi.h>

#include "balanced_hash_relation.h"
#include "balanced_parallel_io.h"
#include "btree/btree_map.h"



enum {STATIC, DYNAMIC};

class parallel_RA
{
private:
    u32 RA_type;
    int* projection_index;

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

    int* get_projection_index()
    {
        return projection_index;
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


    void set_projection_index (s32* index)
    {
        this->projection_index = index;
    }

    u32 get_RA_type() {return RA_type;}



    bool fixed_point_check(relation* output)
    {
        bool fixed_point = false;
        int sum = 0;
        int delta_element_count = output->get_delta_element_count();

        MPI_Allreduce(&delta_element_count, &sum, 1, MPI_INT, MPI_SUM, mcomm.get_comm());

        //std::cout << "delta_element_count " << delta_element_count << "sum " << sum << std::endl;

        if(sum == 0)
            fixed_point = true;

        return fixed_point;
    }


    u32 full_count(relation* output)
    {
        int sum = 0;
        int full_element_count = output->get_full_element_count();
        MPI_Allreduce(&full_element_count, &sum, 1, MPI_INT, MPI_SUM, mcomm.get_comm());
        return sum;
    }



    void full_clique_comm(relation* input, relation* output, u64 *total_buffer_size, u64 **recvbuf)
    {
        u32 buckets = mcomm.get_number_of_buckets();
        u32 nprocs = mcomm.get_nprocs();

        vector_buffer *input_buffer = new vector_buffer[buckets];
        int *input_buffer_size = new int[buckets];

        int* input_distinct_sub_bucket_rank_count;
        int** input_distinct_sub_bucket_rank;
        u32* input_bucket_map;

        int* output_distinct_sub_bucket_rank_count;
        int** output_distinct_sub_bucket_rank;
        u32* output_bucket_map;

        u32** meta_buffer_size = new u32*[buckets];
        memset(meta_buffer_size, 0, sizeof(u32*) * buckets);

        *total_buffer_size = 0;
        u32* bucket_offset = new u32[nprocs];

        /* Meta-data exchange */
        for (u32 i = 0; i < buckets; i++)
        {
            u32 input_arity = input->get_arity();
            //u32 output_arity = output->get_arity();
            //assert (input_arity == 2 && output_arity == 2);
            input_buffer[i] = vector_buffer_create_empty();

            u64 val[input_arity];
            google_relation *full = input->get_full();
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

            input_distinct_sub_bucket_rank_count = input->get_distinct_sub_bucket_rank_count();
            input_distinct_sub_bucket_rank = input->get_distinct_sub_bucket_rank();
            input_bucket_map = input->get_bucket_map();

            output_distinct_sub_bucket_rank_count = output->get_distinct_sub_bucket_rank_count();
            output_distinct_sub_bucket_rank = output->get_distinct_sub_bucket_rank();
            output_bucket_map = output->get_bucket_map();

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
                *total_buffer_size = *total_buffer_size + meta_buffer_size[i][r];

            delete[] req1;
            delete[] stat1;
        }


        /* Allocate buffer */
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

        return;
    }



    void delta_clique_comm(relation* input, relation* output, u64 *total_buffer_size, u64 **recvbuf)
    {
        u32 buckets = mcomm.get_number_of_buckets();
        u32 nprocs = mcomm.get_nprocs();

        //std::cout << "[CLIQUE] buckets " << buckets << std::endl;

        vector_buffer *input_buffer = new vector_buffer[buckets];
        int *input_buffer_size = new int[buckets];

        int* input_distinct_sub_bucket_rank_count;
        int** input_distinct_sub_bucket_rank;
        u32* input_bucket_map;

        int* output_distinct_sub_bucket_rank_count;
        int** output_distinct_sub_bucket_rank;
        u32* output_bucket_map;

        u32** meta_buffer_size = new u32*[buckets];
        memset(meta_buffer_size, 0, sizeof(u32*) * buckets);

        *total_buffer_size = 0;
        u32* bucket_offset = new u32[nprocs];

        u32 input_arity = input->get_arity();
        //u32 output_arity = output->get_arity();
        //assert (input_arity == 2 && output_arity == 2);

        google_relation *delta = input->get_delta();

        input_distinct_sub_bucket_rank_count = input->get_distinct_sub_bucket_rank_count();
        input_distinct_sub_bucket_rank = input->get_distinct_sub_bucket_rank();
        input_bucket_map = input->get_bucket_map();

        output_distinct_sub_bucket_rank_count = output->get_distinct_sub_bucket_rank_count();
        output_distinct_sub_bucket_rank = output->get_distinct_sub_bucket_rank();
        output_bucket_map = output->get_bucket_map();


        /* Meta-data exchange */
        for (u32 i = 0; i < buckets; i++)
        {
            //std::cout << "[CLIQUE] input_distinct_sub_bucket_rank_count " << input_distinct_sub_bucket_rank_count[i] << std::endl;
            input_buffer[i] = vector_buffer_create_empty();
            u64 val[input_arity];

            for (auto it = delta[i].begin(); it != delta[i].end(); it++)
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
            //std::cout << "[CLIQUE] BS: " << input_buffer_size[i] << std::endl;

            meta_buffer_size[i] = new u32[input_distinct_sub_bucket_rank_count[i]];
            memset(meta_buffer_size[i], 0, sizeof(u32) * input_distinct_sub_bucket_rank_count[i]);
#if 1
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
                *total_buffer_size = *total_buffer_size + meta_buffer_size[i][r];

            delete[] req1;
            delete[] stat1;
#endif
        }


        //std::cout << "[CLIQUE] Buffer Size " << *total_buffer_size << std::endl;
#if 1
        /* Allocate buffer */
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
#endif

        delete[] meta_buffer_size;
        delete[] input_buffer;
        delete[] input_buffer_size;
        delete[] bucket_offset;

        return;
    }

    // 1 2  // 2 1
    // 2 3  // 3 2

    int local_join_delta(int input0_buffer_size, int input0_arity, u64 *input0_buffer, relation* i1, relation* output, vector_buffer* local_join_output, int* process_size, int* projection_and_rename)
    {
        u32 local_join_duplicates = 0;
        u32 local_join_inserts = 0;
        google_relation tempT;
        google_relation *input1 = i1->get_delta();
        u32 i1_size = i1->get_delta_element_count();
        //u32 input1_arity = i1->get_arity();
        u32 output_arity = output->get_arity();
        u32* output_sub_bucket_count = output->get_sub_bucket_count();
        u32** output_sub_bucket_rank = output->get_sub_bucket_rank();
        int buckets = mcomm.get_number_of_buckets();
        u64 val[2] = {0, 0};

        //assert(input1_arity == 2 && input0_arity == 2 && output_arity == 2);

        for (int k1 = 0; k1 < input0_buffer_size; k1 = k1 + input0_arity)
        {
            for (int i = 0; i < buckets; i++)
            {
                auto itd = input1[i].find(input0_buffer[k1]);
                if ( itd != input1[i].end() )
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

                                if (projection_and_rename[1] == 0 && projection_and_rename[2] == 1)
                                {
                                    val[0] = it2->first;
                                    val[1] = input0_buffer[k1 + 1];
                                    uint64_t bucket_id = hash_function(val[0]) % buckets;
                                    uint64_t sub_bucket_id = hash_function(val[1]) % output_sub_bucket_count[bucket_id];
                                    int index = output_sub_bucket_rank[bucket_id][sub_bucket_id];
                                    vector_buffer_append(&local_join_output[index], (unsigned char *) val, sizeof(u64) * output_arity);
                                    process_size[index] = process_size[index] + output_arity;
                                }
                                else if (projection_and_rename[1] == 1 && projection_and_rename[2] == 0)
                                {
                                    val[1] = it2->first;
                                    val[0] = input0_buffer[k1 + 1];
                                    uint64_t bucket_id = hash_function(val[0]) % buckets;
                                    uint64_t sub_bucket_id = hash_function(val[1]) % output_sub_bucket_count[bucket_id];
                                    int index = output_sub_bucket_rank[bucket_id][sub_bucket_id];
                                    vector_buffer_append(&local_join_output[index], (unsigned char *) val, sizeof(u64) * output_arity);
                                    process_size[index] = process_size[index] + output_arity;
                                }
                            }
                        }
                        else
                        {
                            Map0* k = new Map0();
                            k->insert(std::make_pair(it2->first, 0));
                            tempT[input0_buffer[k1 + 1]] = k;
                            local_join_inserts++;

                            if (projection_and_rename[1] == 0 && projection_and_rename[2] == 1)
                            {
                                val[0] = it2->first;
                                val[1] = input0_buffer[k1 + 1];
                                uint64_t bucket_id = hash_function(val[0]) % buckets;
                                uint64_t sub_bucket_id = hash_function(val[1]) % output_sub_bucket_count[bucket_id];
                                int index = output_sub_bucket_rank[bucket_id][sub_bucket_id];
                                vector_buffer_append(&local_join_output[index], (unsigned char *) val, sizeof(u64) * output_arity);
                                process_size[index] = process_size[index] + output_arity;
                            }
                            else if (projection_and_rename[1] == 1 && projection_and_rename[2] == 0)
                            {
                                val[1] = it2->first;
                                val[0] = input0_buffer[k1 + 1];
                                uint64_t bucket_id = hash_function(val[0]) % buckets;
                                uint64_t sub_bucket_id = hash_function(val[1]) % output_sub_bucket_count[bucket_id];
                                int index = output_sub_bucket_rank[bucket_id][sub_bucket_id];
                                vector_buffer_append(&local_join_output[index], (unsigned char *) val, sizeof(u64) * output_arity);
                                process_size[index] = process_size[index] + output_arity;
                            }
                        }
                    }
                }
            }
        }


        for(google_relation::iterator ix = tempT.begin(); ix != tempT.end(); ix++)
            delete (ix->second);

        delete[]  input0_buffer;

        std::cout << "[Local Join D] " << i1_size << " " << input0_buffer_size/2 << " " << local_join_inserts << std::endl;
        return local_join_inserts;
    }


    void local_copy(int input0_buffer_size, int input0_arity, u64 *input0_buffer, relation* output, vector_buffer* local_join_output, int* process_size)
    {
        u32 output_arity = output->get_arity();
        u32* output_sub_bucket_count = output->get_sub_bucket_count();
        u32** output_sub_bucket_rank = output->get_sub_bucket_rank();
        int buckets = mcomm.get_number_of_buckets();
        u64 val[2] = {0, 0};

        assert(input0_arity == 2 && output_arity == 2);

        for (int k1 = 0; k1 < input0_buffer_size; k1=k1+input0_arity)
        {
            uint64_t bucket_id = hash_function(input0_buffer[k1+1]) % buckets;
            uint64_t sub_bucket_id = hash_function(input0_buffer[k1]) % output_sub_bucket_count[bucket_id];

            int index = output_sub_bucket_rank[bucket_id][sub_bucket_id];
            val[0] = input0_buffer[k1 + 1];
            val[1] = input0_buffer[k1];
            vector_buffer_append(&local_join_output[index], (unsigned char *) val, sizeof(u64) * output_arity);
            process_size[index] = process_size[index] + output_arity;
        }

        std::cout << "[Local Copy] " << input0_buffer_size/2 << std::endl;

        delete[]  input0_buffer;
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

    // 1 2 // 2 1
    // 2 3 // 3 2

    // 1 3

    int local_inserts_in_full(relation* output)
    {
        //u32 bucket_id = 0;
        int buckets = mcomm.get_number_of_buckets();
        int arity = output->get_arity();
        u32* bucket_map = output->get_bucket_map();
        google_relation *delta = output->get_delta();
        //google_relation *full = output->get_full();
        u32 insert_success = 0;
        u32 insert_attempts = 0;

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

                        //bucket_id = hash_function(t[0]) % buckets;

                        insert_attempts++;
                        if (output->insert_in_full(t) == true)
                            insert_success++;

                        //if (insert_tuple(&(full[bucket_id]), t) == true)
                        //    insert_success++;
                    }
                    delete (local_it->second);
                }

                delta[i].clear();
            }
        }
        output->set_delta_element_count(0);

        std::cout << "[Full Inserts] Attemps " << insert_attempts << " Successful inserts " << insert_success << std::endl;

        return insert_success;
    }

    // 1 2 // 2 1
    // 2 3 // 3 2

    u32 local_join_full(int input0_buffer_size, int input0_arity, u64 *input0_buffer, relation* i1, relation* output, vector_buffer* local_join_output, int* process_size, int* projection_and_rename)
    {
        u32 local_join_duplicates = 0;
        u32 local_join_inserts = 0;
        google_relation tempT;
        google_relation *input1 = i1->get_full();
        u32 i1_size = i1->get_full_element_count();
        u32 output_arity = output->get_arity();
        u32* output_sub_bucket_count = output->get_sub_bucket_count();
        u32** output_sub_bucket_rank = output->get_sub_bucket_rank();
        int buckets = mcomm.get_number_of_buckets();
        u64 val[2] = {0, 0};

        for (int k1 = 0; k1 < input0_buffer_size; k1=k1+input0_arity)
        {
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

                                if (projection_and_rename[1] == 0 && projection_and_rename[2] == 1)
                                {
                                    val[0] = it2->first;
                                    val[1] = input0_buffer[k1 + 1];
                                    uint64_t bucket_id = hash_function(val[0]) % buckets;
                                    uint64_t sub_bucket_id = hash_function(val[1]) % output_sub_bucket_count[bucket_id];
                                    int index = output_sub_bucket_rank[bucket_id][sub_bucket_id];
                                    vector_buffer_append(&(local_join_output[index]), (unsigned char *) val, sizeof(u64) * output_arity);
                                    process_size[index] = process_size[index] + output_arity;
                                }
                                else if (projection_and_rename[1] == 1 && projection_and_rename[2] == 0)
                                {
                                    val[1] = it2->first;
                                    val[0] = input0_buffer[k1 + 1];
                                    uint64_t bucket_id = hash_function(val[0]) % buckets;
                                    uint64_t sub_bucket_id = hash_function(val[1]) % output_sub_bucket_count[bucket_id];
                                    int index = output_sub_bucket_rank[bucket_id][sub_bucket_id];
                                    vector_buffer_append(&(local_join_output[index]), (unsigned char *) val, sizeof(u64) * output_arity);
                                    process_size[index] = process_size[index] + output_arity;
                                }

                            }
                        }
                        else
                        {
                            Map0* k = new Map0();
                            k->insert(std::make_pair(it2->first, 0));
                            tempT[input0_buffer[k1 + 1]] = k;
                            local_join_inserts++;


                            if (projection_and_rename[1] == 0 && projection_and_rename[2] == 1)
                            {
                                val[0] = it2->first;
                                val[1] = input0_buffer[k1 + 1];
                                uint64_t bucket_id = hash_function(val[0]) % buckets;
                                uint64_t sub_bucket_id = hash_function(val[1]) % output_sub_bucket_count[bucket_id];
                                int index = output_sub_bucket_rank[bucket_id][sub_bucket_id];
                                vector_buffer_append(&(local_join_output[index]), (unsigned char *) val, sizeof(u64) * output_arity);
                                process_size[index] = process_size[index] + output_arity;
                            }
                            else if (projection_and_rename[1] == 1 && projection_and_rename[2] == 0)
                            {
                                val[1] = it2->first;
                                val[0] = input0_buffer[k1 + 1];
                                uint64_t bucket_id = hash_function(val[0]) % buckets;
                                uint64_t sub_bucket_id = hash_function(val[1]) % output_sub_bucket_count[bucket_id];
                                int index = output_sub_bucket_rank[bucket_id][sub_bucket_id];
                                vector_buffer_append(&(local_join_output[index]), (unsigned char *) val, sizeof(u64) * output_arity);
                                process_size[index] = process_size[index] + output_arity;
                            }
                        }
                    }
                }
            }
        }

        std::cout << "[Local Join F] " << i1_size << " " << input0_buffer_size/2 << " " << local_join_inserts << std::endl;

        google_relation::iterator ix = tempT.begin();
        for(; ix != tempT.end(); ix++)
            delete (ix->second);
        delete[]  input0_buffer;

        return local_join_inserts;
    }

    void all_to_all(vector_buffer* local_join_output, int* process_size, int* outer_hash_buffer_size, u64** outer_hash_data)
    {
        int nprocs = mcomm.get_nprocs();

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


        /* This step prepares for actual data transfer */
        /* Every process sends to every other process the amount of data it is going to send */

        int recv_process_size_buffer[nprocs];
        memset(recv_process_size_buffer, 0, nprocs * sizeof(int));
        MPI_Alltoall(process_size, 1, MPI_INT, recv_process_size_buffer, 1, MPI_INT, mcomm.get_comm());

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
        //delete[] local_join_output;

    }


    u32 local_insert_in_delta(relation* output, u32 outer_hash_buffer_size, u64* hash_buffer)
    {
        u32 successful_insert = 0;
        u32 output_arity = output->get_arity();
        u64 tuple[2];
        //std::cout << "[WEIRD]" << output->get_full_element_count() << std::endl;
        for (u32 i = 0; i < outer_hash_buffer_size; i = i + output_arity)
        {
            tuple[0] = hash_buffer[i];
            tuple[1] = hash_buffer[i + 1];
            //std::cout << hash_buffer[i] << " " << hash_buffer[i + 1] << std::endl;
            if (output->find_in_full(tuple) == false)
            {
                output->insert_in_delta(tuple);
                successful_insert++;
            }
        }

        delete[] hash_buffer;

        std::cout << "[Delta Inserts] Trying to insert " << outer_hash_buffer_size/2 << " successful_insert " << successful_insert << std::endl;
        return successful_insert;
    }
};
#endif
