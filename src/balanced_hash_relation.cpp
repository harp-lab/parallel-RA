#include "balanced_hash_relation.h"



relation::~relation()
{
    delete[] distinct_sub_bucket_rank_count;
    for (u64 b = 0; b < buckets; b++)
        delete[] distinct_sub_bucket_rank[b];
    delete[] distinct_sub_bucket_rank;


    for (u32 i = 0; i < buckets; i++)
    {
        full[i].remove_tuple();
        delta[i].remove_tuple();
        newt[i].remove_tuple();

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
    delete[] sub_bucket_per_bucket_count;
    delete[] sub_bucket_rank;
}


void relation::local_insert_in_delta()
{
    delete[] delta;
    delta = newt;
    delta_element_count = newt_element_count;
    //std::cout << "delta_element_count " << delta_element_count << std::endl;

    memcpy(delta_bucket_element_count, newt_bucket_element_count, buckets * sizeof(u32));
    for (u32 b = 0; b < buckets; b++)
    {
        memcpy(delta_sub_bucket_element_count[b], newt_sub_bucket_element_count[b], sub_bucket_per_bucket_count[b] * sizeof(u32));
        memset(newt_sub_bucket_element_count[b], 0, sub_bucket_per_bucket_count[b] * sizeof(u32));
    }

    newt = new google_relation[buckets];
    newt_element_count = 0;
    memset(newt_bucket_element_count, 0, buckets * sizeof(u32));

    return;
}

void relation::create_newt()
{
    newt = new google_relation[buckets];
}


void relation::print()
{
    if (mcomm.get_rank() == 0)
    {
        vector_buffer *vb_full = new vector_buffer[buckets];
        std::cout << "FULL ";
        for (u32 i=0; i < buckets; i++)
        {
            vb_full[i] = vector_buffer_create_empty();
            std::vector<u64> prefix = {};
            full[i].as_vector_buffer(&(vb_full[i]), prefix);

            std::cout << vb_full[i].size/sizeof(u64) << " arity " << arity << std::endl;
            for (u32 j=0; j < vb_full[i].size/sizeof(u64); j = j + arity)
            {
                for (int k = 0; k < arity; k++)
                {
                    u64 temp;
                    memcpy(&temp, (vb_full[i].buffer) + (j + k)*sizeof(u64), sizeof(u64));
                    std::cout << temp << " ";
                }
                std::cout << std::endl;
            }

            vector_buffer_free(&(vb_full[i]));
        }
        delete[] vb_full;


        vector_buffer *vb_delta = new vector_buffer[buckets];
        std::cout << "DELTA ";
        for (u32 i=0; i < buckets; i++)
        {
            vb_delta[i] = vector_buffer_create_empty();
            std::vector<u64> prefix = {};
            delta[i].as_vector_buffer(&(vb_delta[i]), prefix);

            std::cout << vb_delta[i].size/sizeof(u64) << " arity " << arity << std::endl;
            for (u32 j=0; j < vb_delta[i].size/sizeof(u64); j = j + arity)
            {
                for (int k = 0; k < arity; k++)
                {
                    u64 temp;
                    memcpy(&temp, (vb_delta[i].buffer) + (j + k)*sizeof(u64), sizeof(u64));
                    std::cout << temp << " ";
                }
                std::cout << std::endl;
            }

            vector_buffer_free(&(vb_delta[i]));
        }
        delete[] vb_delta;


        vector_buffer *vb_newt = new vector_buffer[buckets];
        std::cout << "NEWT ";
        for (u32 i=0; i < buckets; i++)
        {
            vb_newt[i] = vector_buffer_create_empty();
            std::vector<u64> prefix = {};
            newt[i].as_vector_buffer(&(vb_newt[i]), prefix);

            std::cout << vb_newt[i].size/sizeof(u64) << " arity " << arity << std::endl;
            for (u32 j=0; j < vb_newt[i].size/sizeof(u64); j = j + arity)
            {
                for (int k = 0; k < arity; k++)
                {
                    u64 temp;
                    memcpy(&temp, (vb_newt[i].buffer) + (j + k)*sizeof(u64), sizeof(u64));
                    std::cout << temp << " ";
                }
                std::cout << std::endl;
            }

            vector_buffer_free(&(vb_newt[i]));
        }
        delete[] vb_newt;
    }

}


void relation::initialize_relation(u32 arity, mpi_comm& mcomm)
{
    this->arity = arity;
    this->mcomm = mcomm;

    buckets = mcomm.get_nprocs();
    default_sub_bucket_per_bucket_count = 1;
    int rank = mcomm.get_rank();
    int nprocs = mcomm.get_nprocs();

    full_element_count = 0;
    delta_element_count = 0;
    delta = new google_relation[buckets];
    full = new google_relation[buckets];
    newt = new google_relation[buckets];

    sub_bucket_per_bucket_count = new u32[buckets];
    for (u64 b = 0; b < buckets; b++)
    {
        delta[b].set_arity(arity);
        full[b].set_arity(arity);
        newt[b].set_arity(arity);
        sub_bucket_per_bucket_count[b] = default_sub_bucket_per_bucket_count;
    }

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
        sub_bucket_rank[b] = new u32[sub_bucket_per_bucket_count[b]];
        std::unordered_set<int> distinct_ranks;
        for (u64 x = 0; x < sub_bucket_per_bucket_count[b]; x++)
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

        full_sub_bucket_element_count[b] = new u32[sub_bucket_per_bucket_count[b]];
        memset(full_sub_bucket_element_count[b], 0, sizeof(u32) * sub_bucket_per_bucket_count[b]);

        delta_sub_bucket_element_count[b] = new u32[sub_bucket_per_bucket_count[b]];
        memset(delta_sub_bucket_element_count[b], 0, sizeof(u32) * sub_bucket_per_bucket_count[b]);

        newt_sub_bucket_element_count[b] = new u32[sub_bucket_per_bucket_count[b]];
        memset(newt_sub_bucket_element_count[b], 0, sizeof(u32) * sub_bucket_per_bucket_count[b]);
    }
}

void relation::initialize_full(u32 buffer_size, u64 col_count, u64* buffer)
{
    u32 counter = 0;
    u64 t[col_count];
    for (u32 i = 0; i < buffer_size; i = i + arity)
    {
        uint64_t bucket_id = hash_function(buffer[i]) % buckets;

        for (u32 a = i; a < i + arity; a++)
            t[a-i] = buffer[a];
        if (full[bucket_id].insert_tuple_from_array(t, col_count) == true)
        {
            full_element_count++;
            full_bucket_element_count[bucket_id]++;
            counter++;
        }
    }
    //std::cout << "Count = " << counter << std::endl;

    return;
}


void relation::initialize_delta (u32 buffer_size, u64 col_count, u64* buffer)
{
    u64 t[col_count];
    for (u32 i = 0; i < buffer_size; i = i + arity)
    {
        uint64_t bucket_id = hash_function(buffer[i]) % buckets;

        for (u32 a = i; a < i + arity; a++)
            t[a-i] = buffer[a];
        if (delta[bucket_id].insert_tuple_from_array(t, col_count) == true)
        {
            delta_element_count++;
            delta_bucket_element_count[bucket_id]++;
        }
    }

    return;
}

void relation::read_from_file(const char *fname, int full_delta)
{
    /* Parallel I/O to read the relations in parallel */
    file_io.parallel_read_input_relation_from_file_to_local_buffer(fname, mcomm.get_comm());
    file_io.buffer_data_to_hash_buffer_col(buckets, sub_bucket_rank, sub_bucket_per_bucket_count);
    file_io.delete_raw_buffers();

    /* Copy data from buffer to relation */
    if (full_delta == DELTA)
        initialize_delta(file_io.get_hash_buffer_size(), file_io.get_col_count(), file_io.get_hash_buffer());

    else if (full_delta == FULL)
        initialize_full(file_io.get_hash_buffer_size(), file_io.get_col_count(), file_io.get_hash_buffer());


    file_io.delete_hash_buffers();
}


void relation::read_from_relation(relation* input, int full_delta)
{
    google_relation* full = input->get_full();
    std::vector<u64> prefix = {};
    vector_buffer vb = vector_buffer_create_empty();
    full[mcomm.get_rank()].as_vector_buffer(&vb, prefix);

    if (full_delta == DELTA)
        initialize_delta(vb.size/sizeof(u64), 2, (u64*)vb.buffer);

    else if (full_delta == FULL)
        initialize_full(vb.size/sizeof(u64), 2, (u64*)vb.buffer);

    vector_buffer_free(&vb);
}

void relation::flush_full()
{
    google_relation* full = this->get_full();
    full[mcomm.get_rank()].remove_tuple();
    full_element_count = 0;
}


bool relation::find_in_full(u64* t)
{
    uint64_t bucket_id = hash_function(t[0]) % buckets;
    //std::cout << "B ID " << bucket_id << std::endl;
    return full[bucket_id].find_tuple_from_array(t, arity);
}

//bool relation::find_in_full(std::vector<u64> t)
//{
//    uint64_t bucket_id = hash_function(t[0]) % buckets;
//    //std::cout << "B ID " << bucket_id << std::endl;
//    return full[bucket_id].find_tuple_from_vector(t);
//}


bool relation::find_in_delta(u64* t)
{
    uint64_t bucket_id = hash_function(t[0]) % buckets;
    return delta[bucket_id].find_tuple_from_array(t, arity);
}

bool relation::find_in_newt(u64* t)
{
    uint64_t bucket_id = hash_function(t[0]) % buckets;
    //std::cout << "Bucket ID " << bucket_id << std::endl;
    return newt[bucket_id].find_tuple_from_array(t, arity);
}


bool relation::insert_in_delta(u64* t)
{
    uint64_t bucket_id = hash_function(t[0]) % buckets;
    u32 sub_bucket_id = hash_function(t[1]) % sub_bucket_per_bucket_count[bucket_id];

    if (delta[bucket_id].insert_tuple_from_array(t, arity) == true)
    {
        delta_element_count++;
        delta_bucket_element_count[bucket_id]++;
        delta_sub_bucket_element_count[bucket_id][sub_bucket_id]++;
        bucket_map[bucket_id] = 1;

        return true;
    }
    return false;
}


bool relation::insert_in_newt(u64* t)
{
    uint64_t bucket_id = hash_function(t[0]) % buckets;
    u32 sub_bucket_id = hash_function(t[1]) % sub_bucket_per_bucket_count[bucket_id];

    if (newt[bucket_id].insert_tuple_from_array(t, arity) == true)
    {
        newt_element_count++;
        newt_bucket_element_count[bucket_id]++;
        newt_sub_bucket_element_count[bucket_id][sub_bucket_id]++;
        bucket_map[bucket_id] = 1;

        return true;
    }
    return false;
}



bool relation::insert_in_full(u64* t)
{
    uint64_t bucket_id = hash_function(t[0]) % buckets;
    u32 sub_bucket_id = hash_function(t[1]) % sub_bucket_per_bucket_count[bucket_id];


    if (full[bucket_id].insert_tuple_from_array(t, arity) == true)
    {
        full_element_count++;
        full_bucket_element_count[bucket_id]++;
        full_sub_bucket_element_count[bucket_id][sub_bucket_id]++;
        bucket_map[bucket_id] = 1;

        return true;
    }
    return false;
}


int relation::insert_delta_in_full()
{
    u32 insert_success = 0;

    vector_buffer *input_buffer = new vector_buffer[buckets];
    for (u32 i = 0; i < buckets; i++)
    {
        input_buffer[i] = vector_buffer_create_empty();
        if (bucket_map[i] == 1)
        {
            std::vector<u64> prefix = {};
            delta[i].as_vector_buffer(&(input_buffer[i]), prefix);
            for (u64 j = 0; j < (&input_buffer[i])->size / sizeof(u64); j=j+arity)
            {
                //u64 t1, t2;
                //memcpy(&t1, (input_buffer[i].buffer) + (j*sizeof(u64)), sizeof(u64));
                //memcpy(&t2, (input_buffer[i].buffer) + ((j+1)*sizeof(u64)), sizeof(u64));
                //std::cout << "RRRRRRRRRRRRRR " << t1 << " " << t2 << std::endl;
                if (insert_in_full ( (u64*)( (input_buffer[i].buffer) + (j*sizeof(u64)) )) == true)
                    insert_success++;
            }

            //TODO
            delta[i].remove_tuple();

            vector_buffer_free(&input_buffer[i]);
        }
    }
    set_delta_element_count(0);
    set_full_inserts_element_count(insert_success);
    delete[] input_buffer;

    //std::cout << "SUCCESSFUL Inserts " << insert_success << std::endl;

    return insert_success;
}

bool relation::fixed_point_check()
{
    bool fixed_point1 = false;
    int sum = 0;
    MPI_Allreduce(&delta_element_count, &sum, 1, MPI_INT, MPI_SUM, mcomm.get_comm());

    bool fixed_point2 = false;
    int sum1 = 0;
    MPI_Allreduce(&full_inserts_element_count, &sum1, 1, MPI_INT, MPI_SUM, mcomm.get_comm());

    if(sum == 0)
        fixed_point1 = true;

    if(sum1 == 0)
        fixed_point2 = true;

    return fixed_point1 & fixed_point2;
}

#if 0
int relation::load_balance_full()
{

}
#endif
