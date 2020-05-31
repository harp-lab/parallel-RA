#include "parallel_RA_inc.h"


void relation::finalize_relation()
{
    u32 buckets = get_bucket_count();
    newt_element_count = 0;
    full_element_count = 0;
    delta_element_count = 0;

    initailization_type = -1;

    delete[] distinct_sub_bucket_rank_count;
    for (u64 b = 0; b < buckets; b++)
        delete[] distinct_sub_bucket_rank[b];
    delete[] distinct_sub_bucket_rank;

    //std::cout << "Delete bucket count " << buckets << std::endl;

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
    u32 buckets = get_bucket_count();
    delete[] delta;
    delta = newt;
    delta_element_count = newt_element_count;

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

void relation::create_newt()    {    newt = new google_relation[get_bucket_count()];}


void relation::print()
{
    u32 buckets = get_bucket_count();
    if (mcomm.get_rank() == 0)
    {
        vector_buffer *vb_full = new vector_buffer[buckets];
        std::cout << "FULL ";
        for (u32 i=0; i < buckets; i++)
        {
            vb_full[i].vector_buffer_create_empty();
            std::vector<u64> prefix = {};
            full[i].as_vector_buffer_recursive(&(vb_full[i]), prefix);

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

            vb_full[i].vector_buffer_free();
        }
        delete[] vb_full;


        vector_buffer *vb_delta = new vector_buffer[buckets];
        std::cout << "DELTA ";
        for (u32 i=0; i < buckets; i++)
        {
            vb_delta[i].vector_buffer_create_empty();
            std::vector<u64> prefix = {};
            delta[i].as_vector_buffer_recursive(&(vb_delta[i]), prefix);

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

            vb_delta[i].vector_buffer_free();
        }
        delete[] vb_delta;


        vector_buffer *vb_newt = new vector_buffer[buckets];
        std::cout << "NEWT ";
        for (u32 i=0; i < buckets; i++)
        {
            vb_newt[i].vector_buffer_create_empty();
            std::vector<u64> prefix = {};
            newt[i].as_vector_buffer_recursive(&(vb_newt[i]), prefix);

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

            vb_newt[i].vector_buffer_free();
        }
        delete[] vb_newt;
    }

}



void relation::copy_relation(relation*& recv_rel, mpi_comm output_comm, int target_cumulative_rank, int tuples_per_task, u32 input_buckets, u32 output_buckets)
{

#if 1
    /// Meta data Exchange
    //u32 input_buckets = get_bucket_count();
    //u32 output_buckets = output_comm.get_local_nprocs();
    //std::cout << mcomm.get_rank() << " " << target_cumulative_rank << " Input buckets : Output buckets :: " << input_buckets << " : " << output_buckets << std::endl;

    u32 *output_sub_bucket_per_bucket_count = new u32[output_buckets];
    u32 **output_sub_bucket_rank = new u32*[output_buckets];

    for (u32 b = 0; b < output_buckets; b++)
        output_sub_bucket_per_bucket_count[b] = default_sub_bucket_per_bucket_count;

    int rcount = 0;
    for (u32 b = 0; b < output_buckets; b++)
    {
        output_sub_bucket_rank[b] = new u32[output_sub_bucket_per_bucket_count[b]];
        for (u64 x = 0; x < output_sub_bucket_per_bucket_count[b]; x++)
        {
            output_sub_bucket_rank[b][x] = rcount % output_buckets;
            rcount++;
        }
    }

    //google_relation* full = input->get_full();
    vector_buffer *full_input_buffer = new vector_buffer[input_buckets];
    for (u32 i = 0; i < input_buckets; i++)
        full_input_buffer[i].vector_buffer_create_empty();
    int full_process_size[output_comm.get_nprocs()];
    memset(full_process_size, 0, output_comm.get_nprocs() * sizeof(int));

    vector_buffer* full_output = new vector_buffer[output_comm.get_nprocs()];
    for (int i = 0; i < output_comm.get_nprocs(); i++)
        full_output[i].vector_buffer_create_empty();

    //google_relation* delta = input->get_delta();
    vector_buffer *delta_input_buffer = new vector_buffer[input_buckets];
    for (u32 i = 0; i < input_buckets; i++)
        full_input_buffer[i].vector_buffer_create_empty();

    int delta_process_size[output_comm.get_nprocs()];
    memset(delta_process_size, 0, output_comm.get_nprocs() * sizeof(int));

    vector_buffer* delta_output = new vector_buffer[output_comm.get_nprocs()];
    for (int i = 0; i < output_comm.get_nprocs(); i++)
        delta_output[i].vector_buffer_create_empty();


    //MPI_Barrier(MPI_COMM_WORLD);
    //if (mcomm.get_rank() == 0)
    //    std::cout << "A1" << std::endl;
    int fsize = 0;
    int dsize = 0;
    if (tuples_per_task != 0)
    {
        for (u32 i = 0; i < input_buckets; i++)
        {
            std::vector<u64> prefix = {};
            full[i].as_vector_buffer_recursive(&(full_input_buffer[i]), prefix);

            fsize = fsize + (full_input_buffer[i].size / sizeof(u64));
            for (u32 s = 0; s < full_input_buffer[i].size / sizeof(u64); s=s+arity)
            {
                u64 reordered_cur_path[arity];
                for (int j =0; j < arity; j++)
                    memcpy(reordered_cur_path + j, (&full_input_buffer[i])->buffer + ((s + j) * sizeof(u64)), sizeof(u64));

                uint64_t bucket_id = hash_function(reordered_cur_path[0]) % output_buckets;
                uint64_t sub_bucket_id = hash_function(reordered_cur_path[1]) % output_sub_bucket_per_bucket_count[bucket_id];
                int index = output_sub_bucket_rank[bucket_id][sub_bucket_id] + target_cumulative_rank;

                //if (reordered_cur_path[0] == 49 && reordered_cur_path[1] == 47)
                //    std::cout << "YYYYYYYYY " << mcomm.get_rank() << " " << index << " " << target_cumulative_rank << " [" << reordered_cur_path[0] << ", " << reordered_cur_path[1] << "] Bucket ID " << bucket_id << " subbucket id " << sub_bucket_id << " TCR " << target_cumulative_rank << " output_buckets " << output_buckets << " input_buckets " << input_buckets <<  std::endl;

                //if (mcomm.get_rank() == 1)
                //    std::cout << "Index is " << index << " output buckets " << output_buckets << std::endl;
                full_process_size[index] = full_process_size[index] + arity;
                full_output[index].vector_buffer_append((const unsigned char*)reordered_cur_path, sizeof(u64)*arity);
            }
            full_input_buffer[i].vector_buffer_free();
        }

        for (u32 i = 0; i < input_buckets; i++)
        {
            std::vector<u64> prefix = {};
            delta[i].as_vector_buffer_recursive(&(delta_input_buffer[i]), prefix);

            dsize = dsize + (delta_input_buffer[i].size / sizeof(u64));
            for (u32 s = 0; s < delta_input_buffer[i].size / sizeof(u64); s=s+arity)
            {
                u64 reordered_cur_path[arity];
                for (int j =0; j < arity; j++)
                    memcpy(reordered_cur_path + j, (&delta_input_buffer[i])->buffer + ((s + j) * sizeof(u64)), sizeof(u64));

                uint64_t bucket_id = hash_function(reordered_cur_path[0]) % output_buckets;
                uint64_t sub_bucket_id = hash_function(reordered_cur_path[1]) % output_sub_bucket_per_bucket_count[bucket_id];
                int index = output_sub_bucket_rank[bucket_id][sub_bucket_id] + target_cumulative_rank;

                delta_process_size[index] = delta_process_size[index] + arity;
                delta_output[index].vector_buffer_append((const unsigned char*)reordered_cur_path, sizeof(u64)*arity);
            }
            delta_input_buffer[i].vector_buffer_free();
        }
    }

    for (u32 b = 0; b < output_buckets; b++)
        delete[] output_sub_bucket_rank[b];
    delete[] output_sub_bucket_rank;

    delete[] output_sub_bucket_per_bucket_count;


    delete[] delta_input_buffer;
    delete[] full_input_buffer;

    //MPI_Barrier(MPI_COMM_WORLD);
    //if (mcomm.get_rank() == 1)
    //    std::cout << "A2" << std::endl;


#if 1
    u64 full_buffer_size;
    u64* full_buffer;
    all_to_all_comm(full_output, full_process_size, &full_buffer_size, &full_buffer, mcomm.get_comm());
    delete[] full_output;


    finalize_relation();
    recv_rel->set_initailization_type(-1);
    recv_rel->initialize_relation(output_comm);


    u64 t[arity];
    for (u32 i=0; i < full_buffer_size; i=i+arity)
    {
        t[0] = full_buffer[i];
        t[1] = full_buffer[i+1];
        recv_rel->insert_in_full(t);
    }


    u64 delta_buffer_size;
    u64* delta_buffer;
    all_to_all_comm(delta_output, delta_process_size, &delta_buffer_size, &delta_buffer, mcomm.get_comm());
    delete[] delta_output;



    int c = 0;
    for (u32 i=0; i < delta_buffer_size; i=i+arity)
    {
        t[0] = delta_buffer[i];
        t[1] = delta_buffer[i+1];
        if (recv_rel->insert_in_delta(t) == true)
            c++;
    }

    delete[] full_buffer;
    delete[] delta_buffer;

#endif
#endif


    return;
}



void relation::initialize_relation(mpi_comm& mcomm)
{
    this->mcomm = mcomm;

    u32 buckets = mcomm.get_local_nprocs();

    default_sub_bucket_per_bucket_count = 1;
    int rank = mcomm.get_local_rank();
    int nprocs = mcomm.get_local_nprocs();

    //std::cout << "Create Buckets " << buckets << "[" << rank << ", " << nprocs << "]" << std::endl;

    newt_element_count = 0;
    full_element_count = 0;
    delta_element_count = 0;
    delta = new google_relation[buckets];
    full = new google_relation[buckets];
    newt = new google_relation[buckets];

    sub_bucket_per_bucket_count = new u32[buckets];
    for (u32 b = 0; b < buckets; b++)
        sub_bucket_per_bucket_count[b] = default_sub_bucket_per_bucket_count;

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
    for (u32 b = 0; b < buckets; b++)
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


    if (initailization_type != -1)
    {
        //std::cout << "Read from file" << std::endl;
        read_from_file();
    }
}



void relation::initialize(u32 arity)
{
    this->arity = arity;
    full_element_count = 0;
    delta_element_count = 0;
}



void relation::initialize(u32 arity, int tuple_count, char* filename, int version)
{
    if (version == FULL)
    {
        full_element_count = tuple_count;
        delta_element_count = 0;
    }
    else if (version == DELTA)
    {
        delta_element_count = tuple_count;
        full_element_count = 0;
    }

    this->arity = arity;
    this->initailization_type = version;
    this->filename = filename;
}


void relation::initialize_full(u32 buffer_size, u64 col_count, u64* buffer)
{
    u32 counter = 0;
    u64 t[col_count];
    u32 buckets = get_bucket_count();

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
    u32 buckets = get_bucket_count();

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

void relation::read_from_file()
{
    /* Parallel I/O to read the relations in parallel */
    file_io.parallel_read_input_relation_from_file_to_local_buffer(filename, mcomm.get_local_comm());

    file_io.buffer_data_to_hash_buffer_col(get_bucket_count(), sub_bucket_rank, sub_bucket_per_bucket_count, mcomm.get_local_comm());

    file_io.delete_raw_buffers();
#if 1
    /* Copy data from buffer to relation */
    if (initailization_type == DELTA)
        initialize_delta(file_io.get_hash_buffer_size(), file_io.get_col_count(), file_io.get_hash_buffer());

    else if (initailization_type == FULL)
        initialize_full(file_io.get_hash_buffer_size(), file_io.get_col_count(), file_io.get_hash_buffer());

    file_io.delete_hash_buffers();
#endif
}


void relation::read_from_relation(relation* input, int full_delta)
{
    google_relation* full = input->get_full();
    std::vector<u64> prefix = {};
    vector_buffer vb;
    vb.vector_buffer_create_empty();

    full[mcomm.get_rank()].as_vector_buffer_recursive(&vb, prefix);

    if (full_delta == DELTA)
        initialize_delta(vb.size/sizeof(u64), 2, (u64*)vb.buffer);

    else if (full_delta == FULL)
        initialize_full(vb.size/sizeof(u64), 2, (u64*)vb.buffer);

    vb.vector_buffer_free();
}

void relation::flush_full()
{
    google_relation* full = this->get_full();
    full[mcomm.get_rank()].remove_tuple();
    full_element_count = 0;
}


bool relation::find_in_full(u64* t)
{
    uint64_t bucket_id = hash_function(t[0]) % get_bucket_count();
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
    uint64_t bucket_id = hash_function(t[0]) % get_bucket_count();
    return delta[bucket_id].find_tuple_from_array(t, arity);
}

bool relation::find_in_newt(u64* t)
{
    uint64_t bucket_id = hash_function(t[0]) % get_bucket_count();
    //std::cout << "Bucket ID " << bucket_id << std::endl;
    return newt[bucket_id].find_tuple_from_array(t, arity);
}


bool relation::insert_in_delta(u64* t)
{
    uint64_t bucket_id = hash_function(t[0]) % get_bucket_count();
    u32 sub_bucket_id = hash_function(t[1]) % sub_bucket_per_bucket_count[bucket_id];

    //assert((int)bucket_id == mcomm.get_local_rank());
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
    uint64_t bucket_id = hash_function(t[0]) % get_bucket_count();
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
    u32 buckets = get_bucket_count();
    uint64_t bucket_id = hash_function(t[0]) % buckets;
    u32 sub_bucket_id = hash_function(t[1]) % sub_bucket_per_bucket_count[bucket_id];

    //if ((int)bucket_id != mcomm.get_local_rank())
    //    std::cout << "XXXXX [" << t[0] << ", " << t[1] << "] BID " << bucket_id << " [" << mcomm.get_local_rank() << " " << mcomm.get_local_nprocs() << "] [" << mcomm.get_rank() << " " << mcomm.get_nprocs() << "] " << buckets <<  std::endl;
    //assert((int)bucket_id == mcomm.get_local_rank());
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
    u32 buckets = get_bucket_count();

    vector_buffer *input_buffer = new vector_buffer[buckets];
    for (u32 i = 0; i < buckets; i++)
        input_buffer[i].vector_buffer_create_empty();

    for (u32 i = 0; i < buckets; i++)
    {
        if (bucket_map[i] == 1)
        {
            std::vector<u64> prefix = {};
            delta[i].as_vector_buffer_recursive(&(input_buffer[i]), prefix);
            for (u64 j = 0; j < (&input_buffer[i])->size / sizeof(u64); j=j+arity)
            {
                if (insert_in_full ( (u64*)( (input_buffer[i].buffer) + (j*sizeof(u64)) )) == true)
                    insert_success++;
            }

            //TODO
            delta[i].remove_tuple();

            input_buffer[i].vector_buffer_free();
        }
    }
    set_delta_element_count(0);
    delete[] input_buffer;

    //std::cout << "SIZEs " << newt_element_count << " " << delta_element_count << " " << full_element_count << std::endl;

    return insert_success;
}

void relation::fixed_point_check(std::vector<u64>& history)
{
    int sum1 = 0;
    MPI_Allreduce(&delta_element_count, &sum1, 1, MPI_INT, MPI_SUM, mcomm.get_local_comm());

    int sum2 = 0;
    MPI_Allreduce(&full_element_count, &sum2, 1, MPI_INT, MPI_SUM, mcomm.get_local_comm());

    history.push_back(sum1);
    history.push_back(sum2);

    return;
}

#if 0
int relation::load_balance_full()
{

}
#endif
