#include "load_balancer.h"

#if 0
bool load_balancer::load_balance_merge_full(relation* rel, float rf)
{
    u32 rank = (u32)mcomm.get_rank();
    int nprocs = mcomm.get_nprocs();
    int buckets = rel->get_number_of_buckets();

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

    if (mcount > 0.8* buckets)
    {
        if (rank == 0)
            std::cout << "[YES] Bucket Consolidation [" << mcount << " " << 0.8 * buckets << "] Old Sub Bucket count " << total_old_sub_bucket_count << std::endl;
    }
    else
    {
        if (rank == 0)
            std::cout << "[NO] Bucket Consolidation [" << mcount << " " << 0.8 * buckets << "] Old Sub Bucket count " << total_old_sub_bucket_count << std::endl;
        return false;
    }


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

                    total_sub_bucket_size = total_sub_bucket_size + full_sub_bucket_size[i][j];
                }
            }
        }
    }
    int *global_max = new int[buckets];
    memset(global_max, 0, buckets * sizeof(int));
    int average_sub_bucket_size;
    MPI_Allreduce(max_sub_bucket_size, global_max, buckets, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&total_sub_bucket_size, &global_total_sub_bucket_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    delete[] max_sub_bucket_size;
    average_sub_bucket_size = global_total_sub_bucket_size / total_sub_bucket_count;


    //int rcount =  rel->get_last_rank();// sub_bucket_rank[buckets-1][sub_bucket_count[buckets-1] - 1] + 1;
    for (int b = 0; b < buckets; b++)
    {
        if (global_max[b] > average_sub_bucket_size)
            continue;

        global_new_sub_bucket[b] = global_new_sub_bucket[b] / rf;
        if (global_new_sub_bucket[b] == 0)
            global_new_sub_bucket[b] = 1;
        new_sub_bucket_count = new_sub_bucket_count + global_new_sub_bucket[b];

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
                //rcount++;
                //sub_bucket_rank[b][x] = rcount % nprocs;
                //rel->set_last_rank(rcount);

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
    delete[] global_max;


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

bool load_balancer::load_balance_split_full(relation* rel, float rf)
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
            std::cout << "[NO] G RF " << rf << " Bucket Split -- Global Min " << global_min << " Average bucket size " << average_sub_bucket_size << " Global Max [" << global_global_max << " " << average_global_max/buckets  << "] OLD " << old_total_sub_buckets << " NEW " << new_total_sub_buckets << std::endl;
        return false;
    }
    else if (count != 0)
    {
        if (rank == 0)
            std::cout << "[YES] G RF " << rf << " Bucket Split -- Global Min " << global_min << " Average bucket size " << average_sub_bucket_size << " Global Max [" << global_global_max << " " << average_global_max/buckets  << "] OLD " << old_total_sub_buckets << " NEW " << new_total_sub_buckets << std::endl;
    }


    //int rcount = sub_bucket_rank[buckets-1][sub_bucket_count[buckets-1] - 1] + 1;
    int rcount =  rel->get_last_rank();// sub_bucket_rank[buckets-1][sub_bucket_count[buckets-1] - 1] + 1;

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
                    rcount++;
                    sub_bucket_rank[b][x] = rcount % nprocs;
                    rel->set_last_rank(rcount);
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


bool load_balancer::load_balance_merge_full_and_delta(relation* rel, float rf)
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


    if (mcount > 0.8* buckets)
    {
        if (rank == 0)
            std::cout << "[YES] Bucket Consolidation [" << mcount << " " << 0.8 * buckets << "] Old Sub Bucket count " << total_old_sub_bucket_count << std::endl;
    }
    else
    {
        if (rank == 0)
            std::cout << "[NO] Bucket Consolidation [" << mcount << " " << 0.8 * buckets << "] Old Sub Bucket count " << total_old_sub_bucket_count << std::endl;
        return false;
    }


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

                    total_sub_bucket_size = total_sub_bucket_size + full_sub_bucket_size[i][j];
                }
            }
        }
    }
    int *global_max = new int[buckets];
    memset(global_max, 0, buckets * sizeof(int));
    int average_sub_bucket_size;
    MPI_Allreduce(max_sub_bucket_size, global_max, buckets, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&total_sub_bucket_size, &global_total_sub_bucket_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    delete[] max_sub_bucket_size;
    average_sub_bucket_size = global_total_sub_bucket_size / total_sub_bucket_count;


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

    //int rcount =  rel->get_last_rank();// sub_bucket_rank[buckets-1][sub_bucket_count[buckets-1] - 1] + 1;
    for (int b = 0; b < buckets; b++)
    {
        if (global_max[b] > average_sub_bucket_size)
            continue;

        global_new_sub_bucket[b] = global_new_sub_bucket[b] / rf;
        if (global_new_sub_bucket[b] == 0)
            global_new_sub_bucket[b] = 1;
        new_sub_bucket_count = new_sub_bucket_count + global_new_sub_bucket[b];

        if (sub_bucket_count[b] > global_new_sub_bucket[b])
        {
            //check if the max sub bucket is twice the average sub bucket then skip this
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
                //rcount++;
                //sub_bucket_rank[b][x] = rcount % nprocs;
                //rel->set_last_rank(rcount);

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
    delete[] global_max;

    //if (rank == 0)
    //{
    //    for (int b = 0; b < buckets; b++)
    //        for (u64 x = 0; x < global_new_sub_bucket[b]; x++)
    //            std::cout << "(" << b << ", " << x << ") " << sub_bucket_rank[b][x] << "\t";

    //    std::cout << "\n";
    //}

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


bool load_balancer::load_balance_split_full_and_delta(relation* rel, float rf, int rc)
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
    //int min_sub_bucket_size = INT_MAX;
    int *max_sub_bucket_size = new int[buckets];
    memset(max_sub_bucket_size, 0, buckets * sizeof(int));

    int *min_sub_bucket_size = new int[buckets];
    memset(min_sub_bucket_size, 0, buckets * sizeof(int));

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

    if (rc == 0)
    {
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

                        if ((int)full_sub_bucket_size[i][j] < min_sub_bucket_size[i])
                            min_sub_bucket_size[i] = full_sub_bucket_size[i][j];

                        total_sub_bucket_size = total_sub_bucket_size + full_sub_bucket_size[i][j];
                    }
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < buckets; i++)
        {
            total_sub_bucket_count = total_sub_bucket_count + sub_bucket_count[i];
            if (bucket_map[i] == 1)
            {
                for (u32 j = 0; j < sub_bucket_count[i]; j++)
                {
                    if (delta_sub_bucket_size[i][j] != 0)
                    {
                        if ((int)delta_sub_bucket_size[i][j] > max_sub_bucket_size[i])
                            max_sub_bucket_size[i] = delta_sub_bucket_size[i][j];

                        if ((int)full_sub_bucket_size[i][j] < min_sub_bucket_size[i])
                            min_sub_bucket_size[i] = full_sub_bucket_size[i][j];

                        total_sub_bucket_size = total_sub_bucket_size + delta_sub_bucket_size[i][j];
                    }
                }
            }
        }

    }

    int *global_max = new int[buckets];
    memset(global_max, 0, buckets * sizeof(int));

    int *global_min = new int[buckets];
    memset(global_min, 0, buckets * sizeof(int));

    int average_sub_bucket_size;
    MPI_Allreduce(max_sub_bucket_size, global_max, buckets, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(min_sub_bucket_size, global_min, buckets, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&total_sub_bucket_size, &global_total_sub_bucket_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    delete[] max_sub_bucket_size;
    delete[] min_sub_bucket_size;

    average_sub_bucket_size = global_total_sub_bucket_size / total_sub_bucket_count;

    u32 global_new_sub_bucket[buckets];
    memcpy(global_new_sub_bucket, sub_bucket_count, buckets * sizeof(u32));

    int count = 0;
    int old_total_sub_buckets = 0;
    int new_total_sub_buckets = 0;
    int global_global_max = 0;
    int global_global_min = 0;
    int average_global_max = 0;

    for (int i = 0; i < buckets; i++)
    {
        if (global_new_sub_bucket[i] >= (u32)nprocs)
            continue;

        if (global_max[i] > average_sub_bucket_size * rf * 0.8)
        {
            //global_new_sub_bucket[i] = ((ceil((double) global_max[i]/average_sub_bucket_size)) / rf) * rf;
            global_new_sub_bucket[i] = global_new_sub_bucket[i] * rf;
            count++;
        }

        average_global_max = average_global_max + global_max[i];

        if (global_global_max < global_max[i])
            global_global_max = global_max[i];

        if (global_global_min > global_min[i])
            global_global_min = global_min[i];


        old_total_sub_buckets = old_total_sub_buckets + sub_bucket_count[i];
        new_total_sub_buckets = new_total_sub_buckets + global_new_sub_bucket[i];
    }
    delete[] global_max;
    delete[] global_min;

    if (count == 0)
    {
        if (rank == 0)
            std::cout << "[NO] T RF " << rf << " Bucket Split -- " << " Average bucket size " << average_sub_bucket_size << " Global Max Min [" << global_global_max << " " << global_global_min << " " << average_global_max/buckets  << "] OLD " << old_total_sub_buckets << " NEW " << new_total_sub_buckets << std::endl;
        return false;
    }
    else if (count != 0)
    {
        if (rank == 0)
            std::cout << "[YES] T RF " << rf << " Bucket Split -- " << " Average bucket size " << average_sub_bucket_size << " Global Max Min [" << global_global_max << " " << global_global_min << " " << average_global_max/buckets  << "] OLD " << old_total_sub_buckets << " NEW " << new_total_sub_buckets << std::endl;
    }

    int rcount =  rel->get_last_rank();// sub_bucket_rank[buckets-1][sub_bucket_count[buckets-1] - 1] + 1;

    for (int b = 0; b < buckets; b++)
    {
        if (sub_bucket_count[b] < global_new_sub_bucket[b])
        {
            bucket_map[b] = 0;

            std::unordered_set<int> distinct_t_ranks;
            u32* temp = new u32[sub_bucket_count[b]];
            memcpy(temp, sub_bucket_rank[b], sizeof(u32) * sub_bucket_count[b]);

            delete[] sub_bucket_rank[b];
            sub_bucket_rank[b] = new u32[global_new_sub_bucket[b]];

            memcpy(sub_bucket_rank[b], temp, sizeof(int) * sub_bucket_count[b]);
            delete[] temp;

            for (u64 x = 0; x < global_new_sub_bucket[b]; x++)
            {
                if (x >= sub_bucket_count[b])
                {
                    rcount++;
                    sub_bucket_rank[b][x] = rcount % nprocs;
                    rel->set_last_rank(rcount);
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

    //if (rank == 0)
    //{
    //    for (int b = 0; b < buckets; b++)
    //        for (u64 x = 0; x < global_new_sub_bucket[b]; x++)
    //            std::cout << "(" << b << ", " << x << ") " << sub_bucket_rank[b][x] << "\t";

    //    std::cout << "\n";
    //}



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
#endif
