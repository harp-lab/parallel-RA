// Class to manage basic MPI stuff

#ifndef balanced_hash_relation_H
#define balanced_hash_relation_H

enum {DELTA=0, FULL, FULL_AND_DELTA};

enum {COPY=0, JOIN};

#include <iostream>
#include <string>
#include <unordered_set>
#include <mpi.h>

#include "balanced_parallel_io.h"
#include "comm.h"
#include "btree/btree_map.h"


typedef btree::btree_map<u64, u64> Map0;
typedef btree::btree_map<u64, Map0* > google_relation;

bool insert_tuple(google_relation *full, u64* t)
{
    auto it = full->find(t[0]);
    if( it != full->end() )
    {
        auto it2 = (it->second)->find(t[1]);
        if( it2 != (it->second)->end() ) {
            return false;
        }
        else
        {
            (it->second)->insert(std::make_pair(t[1], 0));
            full->insert(std::make_pair(t[0],it->second));
            return true;
        }
    }
    else {
        Map0 *k = new Map0;
        k->insert(std::make_pair(t[1], 0));
        full->insert(std::make_pair(t[0],k));
        return true;
    }
}


class relation
{

private:

    int load_balance_factor;

    bool is_static;
    int join_column_count;
    int arity;

    u32 full_inserts_element_count;

    google_relation *newt;
    u32 newt_element_count;
    u32 **newt_sub_bucket_element_count;
    u32 *newt_bucket_element_count;

    google_relation *full;
    u32 full_element_count;
    u32 **full_sub_bucket_element_count;
    u32 *full_bucket_element_count;   // TODO (implement this carefully)

    google_relation *delta;
    u32 delta_element_count;
    u32 **delta_sub_bucket_element_count;
    u32 *delta_bucket_element_count;   // TODO (implement this carefully)

    u32 *total_sub_bucket_count;
    u32 *sub_bucket_count;          // sub_bucket_count[i] holds the total number of sub-buckets at bucket index i
    u32** sub_bucket_rank;
    int** distinct_sub_bucket_rank;
    int* distinct_sub_bucket_rank_count;

    u32 *bucket_map;

    mpi_comm mcomm;
    parallel_io file_io;

public:
    ~relation()
    {
        destroy();
    }

    u32* get_bucket_map()
    {
        return bucket_map;
    }

    int* get_distinct_sub_bucket_rank_count()
    {
        return distinct_sub_bucket_rank_count;
    }

    int** get_distinct_sub_bucket_rank()
    {
        return distinct_sub_bucket_rank;
    }


    u32** get_full_sub_bucket_element_count()
    {
        return full_sub_bucket_element_count;
    }


    u32** get_delta_sub_bucket_element_count()
    {
        return delta_sub_bucket_element_count;
    }

    u32** get_new_sub_bucket_element_count()
    {
        return newt_sub_bucket_element_count;
    }


    u32* get_sub_bucket_count()
    {
        return sub_bucket_count;
    }

    u32** get_sub_bucket_rank()
    {
        return sub_bucket_rank;
    }

    u32 get_arity ()
    {
        return arity;
    }


    google_relation* get_full()
    {
        return full;
    }


    google_relation* get_newt()
    {
        return newt;
    }


    void set_delta_element_count(int val)
    {
        //std::cout << "CALLED " << val << std::endl;
        delta_element_count = val;
    }

    int get_delta_element_count()
    {
        return delta_element_count;
    }


    int get_new_element_count()
    {
        return newt_element_count;
    }


    int get_full_inserts_element_count()
    {
        return full_inserts_element_count;
    }

    void set_full_inserts_element_count(int x)
    {
        full_inserts_element_count = x;
    }


    int get_full_element_count()
    {
        return full_element_count;
    }

    google_relation* get_delta()
    {
        return delta;
    }


    void copy_newt_to_delta()
    {
        delta = newt;
    }

    void local_insert_in_delta()
    {
        int rank = mcomm.get_rank();
        delete[] delta;
        delta = newt;
        delta_element_count = newt_element_count;

        u32 buckets = mcomm.get_number_of_buckets();
        memcpy(delta_bucket_element_count, newt_bucket_element_count, buckets * sizeof(u32));
        for (u32 b = 0; b < buckets; b++)
        {
            memcpy(delta_sub_bucket_element_count[b], newt_sub_bucket_element_count[b], sub_bucket_count[b] * sizeof(u32));
            memset(newt_sub_bucket_element_count[b], 0, sub_bucket_count[b] * sizeof(u32));
        }

        newt = new google_relation[buckets];
        newt_element_count = 0;
        memset(newt_bucket_element_count, 0, buckets * sizeof(u32));

        if (rank == 0)
            std::cout << "[Pointer copy from new to delta] " << delta_element_count << std::endl;

        return;
    }



    void create_newt()
    {
        int buckets = mcomm.get_number_of_buckets();
        newt = new google_relation[buckets];
    }


    void delete_delta()
    {
        delete[] delta;
    }


    int get_nprocs()
    {
        return mcomm.get_nprocs();
    }


    void initialize_empty(u32 arity, u32 join_column_count, bool is_static, u32 sbc, mpi_comm mcomm)
    {
        assert(arity == 2);
        this->arity = arity;
        this->is_static = is_static;
        this->join_column_count = join_column_count;
        this->mcomm = mcomm;

        u32 buckets = mcomm.get_number_of_buckets();
        int rank = mcomm.get_rank();
        int nprocs = mcomm.get_nprocs();

        full_element_count = 0;
        delta_element_count = 0;
        delta = new google_relation[buckets];
        full = new google_relation[buckets];
        newt = new google_relation[buckets];

        sub_bucket_count = new u32[buckets];
        for (u64 b = 0; b < buckets; b++)
            sub_bucket_count[b] = sbc;

        sub_bucket_rank = new u32*[buckets];
        distinct_sub_bucket_rank = new int*[buckets];
        distinct_sub_bucket_rank_count = new int[buckets];

        bucket_map = new u32[buckets];
        memset(bucket_map, 0, sizeof(u32) * buckets);

        full_sub_bucket_element_count = new u32*[buckets];
        memset(full_sub_bucket_element_count, 0, sizeof(u32*) * buckets);

        delta_sub_bucket_element_count = new u32*[buckets];
        memset(delta_sub_bucket_element_count, 0, sizeof(u32*) * buckets);

        newt_sub_bucket_element_count = new u32*[buckets];
        memset(newt_sub_bucket_element_count, 0, sizeof(u32*) * buckets);

        full_bucket_element_count = new u32[buckets];
        memset(full_bucket_element_count, 0, sizeof(u32) * buckets);

        delta_bucket_element_count = new u32[buckets];
        memset(delta_bucket_element_count, 0, sizeof(u32) * buckets);

        newt_bucket_element_count = new u32[buckets];
        memset(newt_bucket_element_count, 0, sizeof(u32) * buckets);

        int rcount = 0;
        for (u64 b = 0; b < buckets; b++)
        {
            sub_bucket_rank[b] = new u32[sub_bucket_count[b]];
            std::unordered_set<int> distinct_ranks;
            for (u64 x = 0; x < sub_bucket_count[b]; x++)
            {
                sub_bucket_rank[b][x] = rcount % nprocs;

                if (sub_bucket_rank[b][x] == (u32)rank)
                    bucket_map[b] = 1;

                distinct_ranks.insert(sub_bucket_rank[b][x]);
                rcount++;
            }

            distinct_sub_bucket_rank_count[b] = distinct_ranks.size();
            distinct_sub_bucket_rank[b] = new int[distinct_sub_bucket_rank_count[b]];
            u32 x = 0;
            for (auto it = distinct_ranks.begin(); it != distinct_ranks.end(); ++it, x++ )
                distinct_sub_bucket_rank[b][x] = *it;

            full_sub_bucket_element_count[b] = new u32[sub_bucket_count[b]];
            memset(full_sub_bucket_element_count[b], 0, sizeof(u32) * sub_bucket_count[b]);

            delta_sub_bucket_element_count[b] = new u32[sub_bucket_count[b]];
            memset(delta_sub_bucket_element_count[b], 0, sizeof(u32) * sub_bucket_count[b]);

            newt_sub_bucket_element_count[b] = new u32[sub_bucket_count[b]];
            memset(newt_sub_bucket_element_count[b], 0, sizeof(u32) * sub_bucket_count[b]);
        }
    }


    void initialize(u32 arity, u32 join_column_count, bool is_static, u32 sbc, const char *fname, int full_delta, mpi_comm& mcomm, int col_index)
    {
        assert(arity == 2);
        this->arity = arity;
        this->is_static = is_static;
        this->join_column_count = join_column_count;
        this->mcomm = mcomm;

        u32 buckets = mcomm.get_number_of_buckets();
        int rank = mcomm.get_rank();
        int nprocs = mcomm.get_nprocs();

        //std::cout << "Initialize " << buckets << std::endl;


        full_element_count = 0;
        delta_element_count = 0;
        delta = new google_relation[buckets];
        full = new google_relation[buckets];
        newt = new google_relation[buckets];
        /*
        for (u32 i = 0; i < buckets; i++)
            full->setArity(arity);
        */

        sub_bucket_count = new u32[buckets];
        for (u64 b = 0; b < buckets; b++)
            sub_bucket_count[b] = sbc;

        sub_bucket_rank = new u32*[buckets];
        distinct_sub_bucket_rank = new int*[buckets];
        distinct_sub_bucket_rank_count = new int[buckets];

        bucket_map = new u32[buckets];
        memset(bucket_map, 0, sizeof(u32) * buckets);

        full_sub_bucket_element_count = new u32*[buckets];
        memset(full_sub_bucket_element_count, 0, sizeof(u32*) * buckets);

        delta_sub_bucket_element_count = new u32*[buckets];
        memset(delta_sub_bucket_element_count, 0, sizeof(u32*) * buckets);

        newt_sub_bucket_element_count = new u32*[buckets];
        memset(newt_sub_bucket_element_count, 0, sizeof(u32*) * buckets);

        full_bucket_element_count = new u32[buckets];
        memset(full_bucket_element_count, 0, sizeof(u32) * buckets);

        delta_bucket_element_count = new u32[buckets];
        memset(delta_bucket_element_count, 0, sizeof(u32) * buckets);

        newt_bucket_element_count = new u32[buckets];
        memset(newt_bucket_element_count, 0, sizeof(u32) * buckets);

        int rcount = 0;
        for (u64 b = 0; b < buckets; b++)
        {
            sub_bucket_rank[b] = new u32[sub_bucket_count[b]];
            std::unordered_set<int> distinct_ranks;
            for (u64 x = 0; x < sub_bucket_count[b]; x++)
            {
                sub_bucket_rank[b][x] = rcount % nprocs;

                if (sub_bucket_rank[b][x] == (u32)rank)
                    bucket_map[b] = 1;

                distinct_ranks.insert(sub_bucket_rank[b][x]);
                rcount++;
            }

            distinct_sub_bucket_rank_count[b] = distinct_ranks.size();
            distinct_sub_bucket_rank[b] = new int[distinct_sub_bucket_rank_count[b]];
            u32 x  = 0;
            for ( auto it = distinct_ranks.begin(); it != distinct_ranks.end(); ++it, x++ )
                distinct_sub_bucket_rank[b][x] = *it;

            full_sub_bucket_element_count[b] = new u32[sub_bucket_count[b]];
            memset(full_sub_bucket_element_count[b], 0, sizeof(u32) * sub_bucket_count[b]);

            delta_sub_bucket_element_count[b] = new u32[sub_bucket_count[b]];
            memset(delta_sub_bucket_element_count[b], 0, sizeof(u32) * sub_bucket_count[b]);

            newt_sub_bucket_element_count[b] = new u32[sub_bucket_count[b]];
            memset(newt_sub_bucket_element_count[b], 0, sizeof(u32) * sub_bucket_count[b]);
        }

        /* Parallel I/O to read the relations in parallel */
        file_io.parallel_read_input_relation_from_file_to_local_buffer(fname, mcomm.get_rank(), mcomm.get_nprocs());
        file_io.buffer_data_to_hash_buffer_col(mcomm.get_nprocs(), mcomm.get_comm(), buckets, sub_bucket_rank, sub_bucket_count, col_index);
        file_io.delete_raw_buffers();

        /* Copy data from buffer to relation */
        if (full_delta == DELTA)
            initialize_delta(file_io.get_hash_buffer_size(), file_io.get_hash_buffer());
        else if (full_delta == FULL)
        {
            //std::cout << "COL Count " << file_io.get_col_count() << " ROW Count " << file_io.get_hash_buffer_size() << std::endl;
            initialize_full(file_io.get_hash_buffer_size(), file_io.get_hash_buffer());
        }
        else if (full_delta == FULL_AND_DELTA)
            initialize_full_and_delta(file_io.get_hash_buffer_size(), file_io.get_hash_buffer());

        /* load balance the relation */
        //load_balance(load_balance_factor);

        file_io.delete_hash_buffers();
    }


    void initialize_with_rename_and_projection(u32 arity, u32 join_column_count, bool is_static, u32 sbc, const char *fname, int full_delta, int* rename_and_project_copy, mpi_comm mcomm, int col_index)
    {
        assert(arity == 2);
        this->arity = arity;
        this->is_static = is_static;
        this->join_column_count = join_column_count;
        this->mcomm = mcomm;

        u32 buckets = mcomm.get_number_of_buckets();
        int rank = mcomm.get_rank();
        int nprocs = mcomm.get_nprocs();

        //std::cout << "buckets buckets" << buckets << "buckets" << std::endl;

        full_element_count = 0;
        delta_element_count = 0;
        delta = new google_relation[buckets];
        full = new google_relation[buckets];
        newt = new google_relation[buckets];

        sub_bucket_count = new u32[buckets];
        for (u64 b = 0; b < buckets; b++)
            sub_bucket_count[b] = sbc;

        sub_bucket_rank = new u32*[buckets];
        distinct_sub_bucket_rank = new int*[buckets];
        distinct_sub_bucket_rank_count = new int[buckets];

        bucket_map = new u32[buckets];
        memset(bucket_map, 0, sizeof(u32) * buckets);

        full_sub_bucket_element_count = new u32*[buckets];
        memset(full_sub_bucket_element_count, 0, sizeof(u32*) * buckets);

        delta_sub_bucket_element_count = new u32*[buckets];
        memset(delta_sub_bucket_element_count, 0, sizeof(u32*) * buckets);

        newt_sub_bucket_element_count = new u32*[buckets];
        memset(newt_sub_bucket_element_count, 0, sizeof(u32*) * buckets);

        full_bucket_element_count = new u32[buckets];
        memset(full_bucket_element_count, 0, sizeof(u32) * buckets);

        delta_bucket_element_count = new u32[buckets];
        memset(delta_bucket_element_count, 0, sizeof(u32) * buckets);

        newt_bucket_element_count = new u32[buckets];
        memset(newt_bucket_element_count, 0, sizeof(u32) * buckets);

        int rcount = 0;
        for (u64 b = 0; b < buckets; b++)
        {
            sub_bucket_rank[b] = new u32[sub_bucket_count[b]];
            std::unordered_set<int> distinct_ranks;
            for (u64 x = 0; x < sub_bucket_count[b]; x++)
            {
                sub_bucket_rank[b][x] = rcount % nprocs;

                if (sub_bucket_rank[b][x] == (u32)rank)
                    bucket_map[b] = 1;

                distinct_ranks.insert(sub_bucket_rank[b][x]);
                rcount++;
            }

            distinct_sub_bucket_rank_count[b] = distinct_ranks.size();
            distinct_sub_bucket_rank[b] = new int[distinct_sub_bucket_rank_count[b]];
            u32 x  = 0;
            for ( auto it = distinct_ranks.begin(); it != distinct_ranks.end(); ++it, x++ )
                distinct_sub_bucket_rank[b][x] = *it;

            full_sub_bucket_element_count[b] = new u32[sub_bucket_count[b]];
            memset(full_sub_bucket_element_count[b], 0, sizeof(u32) * sub_bucket_count[b]);

            delta_sub_bucket_element_count[b] = new u32[sub_bucket_count[b]];
            memset(delta_sub_bucket_element_count[b], 0, sizeof(u32) * sub_bucket_count[b]);

            newt_sub_bucket_element_count[b] = new u32[sub_bucket_count[b]];
            memset(newt_sub_bucket_element_count[b], 0, sizeof(u32) * sub_bucket_count[b]);
        }

        /* Parallel I/O to read the relations in parallel */
        file_io.parallel_read_input_relation_from_file_to_local_buffer(fname, mcomm.get_rank(), mcomm.get_nprocs());
        file_io.buffer_data_to_hash_buffer_col(mcomm.get_nprocs(), mcomm.get_comm(), buckets, sub_bucket_rank, sub_bucket_count, col_index);
        file_io.delete_raw_buffers();

        //std::cout << "[INIIIIIIIII] Process " << rank << " " << file_io.get_hash_buffer_size() << std::endl;

        /* Copy data from buffer to relation */
        if (full_delta == DELTA)
            initialize_delta_with_rename_and_project(file_io.get_hash_buffer_size(), file_io.get_hash_buffer(), rename_and_project_copy);

#if 0
        else if (full_delta == FULL)
            initialize_full_with_rename_and_project(file_io.get_col_count() * file_io.get_row_count(), file_io.get_hash_buffer(), rename_and_project_copy);
        else if (full_delta == FULL_AND_DELTA)
            initialize_full_and_delta_with_rename_and_project(file_io.get_col_count() * file_io.get_row_count(), file_io.get_hash_buffer(), rename_and_project_copy);
#endif
        /* load balance the relation */
        //load_balance(load_balance_factor);

        //std::cout << "XXXXXXXXXXXXXXX " << get_delta_element_count() << std::endl;

        file_io.delete_hash_buffers();
    }




    void destroy()
    {
        u32 buckets = mcomm.get_number_of_buckets();
        //std::cout << "Number of buckets " << buckets << std::endl;


        delete[] distinct_sub_bucket_rank_count;
        for (u64 b = 0; b < buckets; b++)
            delete[] distinct_sub_bucket_rank[b];
        delete[] distinct_sub_bucket_rank;


        for (u32 i = 0; i < buckets; i++)
        {
            for(google_relation::iterator ix = full[i].begin(); ix != full[i].end(); ix++)
                delete (ix->second);

            for(google_relation::iterator ix = delta[i].begin(); ix != delta[i].end(); ix++)
                delete (ix->second);

            for(google_relation::iterator ix = newt[i].begin(); ix != newt[i].end(); ix++)
                delete (ix->second);

            delete[] sub_bucket_rank[i];
            delete[] full_sub_bucket_element_count[i];
            delete[] delta_sub_bucket_element_count[i];
            delete[] newt_sub_bucket_element_count[i];
        }

        delete[] delta_sub_bucket_element_count;
        delete[] delta_bucket_element_count;

        delete[] full_bucket_element_count;
        delete[] full_sub_bucket_element_count;

        delete[] newt_bucket_element_count;
        delete[] newt_sub_bucket_element_count;

        delete[] delta;
        delete[] full;
        delete[] newt;
        delete[] bucket_map;
        delete[] sub_bucket_count;
        delete[] sub_bucket_rank;
    }


    void initialize_full (u32 buffer_size, u64* buffer)
    {
        u64 t[arity];
        assert(arity == 2);
        //std::cout << "Buckets: " << buckets << std::endl;
        //std::cout << "Buffer Size: " << buffer_size << std::endl;

        for (u32 i = 0; i < buffer_size; i = i + arity)
        {
            for (u32 a = i; a < i + arity; a++)
                t[a-i] = buffer[a];

            insert_in_full(t);
        }
        //std::cout << "Finish Looping " << full_element_count << std::endl;
    }


    void initialize_delta (u32 buffer_size, u64* buffer)
    {
        u64 t[arity];
        u32 buckets = mcomm.get_number_of_buckets();

        for (u32 i = 0; i < buffer_size; i = i + arity)
        {
            uint64_t bucket_id = hash_function(buffer[i]) % buckets;

            for (u32 a = i; a < i + arity; a++)
                t[a-i] = buffer[a];

            insert_tuple(&(delta[bucket_id]), t);
        }
    }


    void initialize_full_and_delta (u32 buffer_size, u64* buffer)
    {
        u64 t[arity];
        u32 buckets = mcomm.get_number_of_buckets();

        for (u32 i = 0; i < buffer_size; i = i + arity)
        {
            uint64_t bucket_id = hash_function(buffer[i]) % buckets;
            uint64_t sub_bucket_id = hash_function(buffer[i + 1]) % sub_bucket_count[bucket_id];
            bucket_map[bucket_id] = 1;

            for (u32 a = i; a < i + arity; a++)
                t[a-i] = buffer[a];

            if (insert_tuple(&(full[bucket_id]), t) == true)
            {
                insert_tuple(&(delta[bucket_id]), t);
                full_sub_bucket_element_count[bucket_id][sub_bucket_id]++;
                full_bucket_element_count[bucket_id]++;
            }
        }
    }


    void initialize_delta_with_rename_and_project (u32 buffer_size, u64* buffer, int* rename_and_project_copy)
    {
        u64 t[arity];
        //u32 buckets = mcomm.get_number_of_buckets();
        //assert(rename_and_project_copy[0] == 1 && rename_and_project_copy[1] == 0);
        //assert(arity == 2);

        //int rank = mcomm.get_rank();
        if(rename_and_project_copy[0] == 1 && rename_and_project_copy[1] == 0)
        {
            for (u32 i = 0; i < buffer_size; i = i + arity)
            {
                //uint64_t bucket_id = hash_function(buffer[i + 1]) % buckets;
                //uint64_t sub_bucket_id = hash_function(buffer[i]) % sub_bucket_count[bucket_id];


                t[0] = buffer[i + 1];
                t[1] = buffer[i];

                //if (rank == 0)
                //    std::cout << i << " Inserts: " << t[0] << " " << t[1] << std::endl;
                insert_in_delta(t);

                //delta_element_count++;
                //delta_bucket_element_count[bucket_id]++;
                //delta_sub_bucket_element_count[bucket_id][sub_bucket_id]++;
            }
        }
        else
        {
            for (u32 i = 0; i < buffer_size; i = i + arity)
            {
                //uint64_t bucket_id = hash_function(buffer[i]) % buckets;
                //uint64_t sub_bucket_id = hash_function(buffer[i + 1]) % sub_bucket_count[bucket_id];

                t[0] = buffer[i];
                t[1] = buffer[i + 1];
                insert_in_delta(t) == true;

                //if (insert_tuple(&(delta[bucket_id]), t) == true)
                //{
                //    delta_element_count++;
                //    delta_bucket_element_count[bucket_id]++;
                //    delta_sub_bucket_element_count[bucket_id][sub_bucket_id]++;
                //}
            }
        }

        //std::cout << "[INITIALIZE] Elements inserted in delta " << delta_element_count << std::endl << std::endl;
    }


    bool check_for_fixed_point()
    {
        bool fixed_point = false;
        int sum = 0;
        MPI_Allreduce(&delta_element_count, &sum, 1, MPI_INT, MPI_SUM, mcomm.get_comm());

        //std::cout << "delta_element_count " << delta_element_count << " sum " << sum << std::endl;

        if(sum == 0)
            fixed_point = true;

        return fixed_point;
    }


    bool find_in_delta(u32 bucket_id, u64* t)
    {
        auto it = delta[bucket_id].find(t[0]);
        if( it != delta[bucket_id].end() )
        {
            auto it2 = (it->second)->find(t[1]);
            if( it2 != (it->second)->end() )
                return false;
            else
                return true;
        }
        else
            return true;
    }


    bool find_in_full(u64* t)
    {
        u32 buckets = mcomm.get_number_of_buckets();
        uint64_t bucket_id = hash_function(t[0]) % buckets;
        auto it = full[bucket_id].find(t[0]);
        if( it != full[bucket_id].end() )
        {
            auto it2 = (it->second)->find(t[1]);
            if( it2 != (it->second)->end() )
                return true;
            else
                return false;
        }
        else
            return false;
    }


    bool find_in_delta(u64* t)
    {
        u32 buckets = mcomm.get_number_of_buckets();
        uint64_t bucket_id = hash_function(t[0]) % buckets;
        auto it = delta[bucket_id].find(t[0]);
        if( it != delta[bucket_id].end() )
        {
            auto it2 = (it->second)->find(t[1]);
            if( it2 != (it->second)->end() )
                return false;
            else
                return true;
        }
        else
            return true;

    }


    bool insert_in_delta(u64* t)
    {
        u32 buckets = mcomm.get_number_of_buckets();
        uint64_t bucket_id = hash_function(t[0]) % buckets;
        u32 sub_bucket_id = hash_function(t[1]) % sub_bucket_count[bucket_id];
        auto it = delta[bucket_id].find(t[0]);
        if( it != delta[bucket_id].end() )
        {
            auto it2 = (it->second)->find(t[1]);
            if( it2 != (it->second)->end() ) {
                return false;
            }
            else
            {
                (it->second)->insert(std::make_pair(t[1], 0));
                delta[bucket_id].insert(std::make_pair(t[0],it->second));
                delta_element_count++;
                delta_bucket_element_count[bucket_id]++;
                delta_sub_bucket_element_count[bucket_id][sub_bucket_id]++;
                bucket_map[bucket_id] = 1;
                return true;
            }
        }
        else {
            Map0 *k = new Map0;
            k->insert(std::make_pair(t[1], 0));
            delta[bucket_id].insert(std::make_pair(t[0],k));
            delta_element_count++;
            delta_bucket_element_count[bucket_id]++;
            delta_sub_bucket_element_count[bucket_id][sub_bucket_id]++;
            bucket_map[bucket_id] = 1;
            return true;
        }
    }


    bool insert_in_newt(u64* t)
    {
        u32 buckets = mcomm.get_number_of_buckets();
        uint64_t bucket_id = hash_function(t[0]) % buckets;
        u32 sub_bucket_id = hash_function(t[1]) % sub_bucket_count[bucket_id];

        auto it = newt[bucket_id].find(t[0]);
        if( it != newt[bucket_id].end() )
        {
            auto it2 = (it->second)->find(t[1]);
            if( it2 != (it->second)->end() ) {
                return false;
            }
            else
            {
                (it->second)->insert(std::make_pair(t[1], 0));
                newt[bucket_id].insert(std::make_pair(t[0],it->second));
                newt_element_count++;
                newt_bucket_element_count[bucket_id]++;
                newt_sub_bucket_element_count[bucket_id][sub_bucket_id]++;
                bucket_map[bucket_id] = 1;
                return true;
            }
        }
        else {
            Map0 *k = new Map0;
            k->insert(std::make_pair(t[1], 0));
            newt[bucket_id].insert(std::make_pair(t[0],k));
            newt_element_count++;
            newt_bucket_element_count[bucket_id]++;
            newt_sub_bucket_element_count[bucket_id][sub_bucket_id]++;
            bucket_map[bucket_id] = 1;
            return true;
        }
    }



    bool insert_in_full(u64* t)
    {

        u32 buckets = mcomm.get_number_of_buckets();
        uint64_t bucket_id = hash_function(t[0]) % buckets;
        u32 sub_bucket_id = hash_function(t[1]) % sub_bucket_count[bucket_id];
        auto it = full[bucket_id].find(t[0]);
        if( it != full[bucket_id].end() )
        {
            auto it2 = (it->second)->find(t[1]);
            if( it2 != (it->second)->end() ) {
                return false;
            }
            else
            {
                (it->second)->insert(std::make_pair(t[1], 0));
                full[bucket_id].insert(std::make_pair(t[0],it->second));
                full_element_count++;
                full_bucket_element_count[bucket_id]++;
                full_sub_bucket_element_count[bucket_id][sub_bucket_id]++;
                bucket_map[bucket_id] = 1;
                return true;
            }
        }
        else {
            Map0 *k = new Map0;
            k->insert(std::make_pair(t[1], 0));
            full[bucket_id].insert(std::make_pair(t[0],k));
            full_element_count++;
            full_bucket_element_count[bucket_id]++;
            full_sub_bucket_element_count[bucket_id][sub_bucket_id]++;
            bucket_map[bucket_id] = 1;
            return true;
        }
    }


    int insert_delta_in_full()
    {
        int rank = mcomm.get_rank();
        int buckets = mcomm.get_number_of_buckets();
        u32 insert_success = 0;
        u32 insert_attempts = 0;

        u64 t[arity];
        //int i = rank;
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

                        insert_attempts++;
                        if (insert_in_full(t) == true)
                            insert_success++;
                    }
                    delete (local_it->second);
                }
                delta[i].clear();
            }
        }
        set_delta_element_count(0);
        set_full_inserts_element_count(insert_success);

        if (rank == 0)
            std::cout << "[Full Inserts] Attemps " << insert_attempts << " Successful inserts " << insert_success << std::endl;

        return insert_success;
    }
};

#endif
