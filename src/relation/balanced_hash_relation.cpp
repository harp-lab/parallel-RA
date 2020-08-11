/*
 * Relation class
 * Copyright (c) Sidharth Kumar, et al, see License.md
 */



#include "../parallel_RA_inc.h"

u32 relation::get_global_delta_element_count()
{
    int dec = (int)delta_element_count;
    int global_delta_element_count;
    MPI_Allreduce(&dec, &global_delta_element_count, 1, MPI_INT, MPI_SUM, mcomm.get_local_comm());
    return global_delta_element_count;
}


u32 relation::get_global_full_element_count()
{
    u32 global_full_element_count;
    MPI_Allreduce(&full_element_count, &global_full_element_count, 1, MPI_INT, MPI_SUM, mcomm.get_local_comm());
    return global_full_element_count;
}




void relation::print()
{
    u32 buckets = get_bucket_count();
    //if (mcomm.get_rank() == 0)
    {
        vector_buffer *vb_full = new vector_buffer[buckets];
        //std::cout << "FULL ";
        for (u32 i=0; i < buckets; i++)
        {
            vb_full[i].vector_buffer_create_empty();
            std::vector<u64> prefix = {};
            full[i].as_vector_buffer_recursive(&(vb_full[i]), prefix);

            //std::cout << get_debug_id() << " " << mcomm.get_rank() << " FULL " << vb_full[i].size/(sizeof(u64) * (arity + 1)) << " XXXX arity " << arity + 1 << std::endl;

            for (u32 j=0; j < vb_full[i].size/sizeof(u64); j = j + arity+1)
            {
                if (j % (arity+1) == 0)
                    std::cout << "F [" << j/(arity+1) << "] ";
                for (u32 k = 0; k < arity+1; k++)
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

        u32 i = mcomm.get_rank();
        {
            vb_delta[i].vector_buffer_create_empty();
            std::vector<u64> prefix = {};
            delta[i].as_vector_buffer_recursive(&(vb_delta[i]), prefix);

            std::cout << get_debug_id() << " " << mcomm.get_rank() << " DELTA " << vb_delta[i].size/(sizeof(u64) * (arity + 1)) << " arity " << arity + 1 << std::endl;

            for (u32 j=0; j < vb_delta[i].size/sizeof(u64); j = j + arity+1)
            {
                if (j % (arity+1) == 0)
                    std::cout << "D ";

                for (u32 k = 0; k < arity+1; k++)
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

        /*
        vector_buffer *vb_newt = new vector_buffer[buckets];
        std::cout << "NEWT ";
        for (u32 i=0; i < buckets; i++)
        {
            vb_newt[i].vector_buffer_create_empty();
            std::vector<u64> prefix = {};
            newt[i].as_vector_buffer_recursive(&(vb_newt[i]), prefix);

            std::cout << vb_newt[i].size/sizeof(u64) << " arity " << arity+1 << std::endl;
            for (u32 j=0; j < vb_newt[i].size/sizeof(u64); j = j + arity+1)
            {
                if (j % (arity+1) == 0)
                    std::cout << "N ";

                for (u32 k = 0; k < arity+1; k++)
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
        */
    }
}


#if 0
void relation::flush_full()
{
    full[mcomm.get_rank()]->remove_tuple();
    full_element_count = 0;
    for (int i =0; i < get_bucket_count(); i++)
    {
        full_bucket_element_count[i] = 0;
        for (u32 j=0; j < sub_bucket_per_bucket_count[i]; j++)
            full_sub_bucket_element_count[i][j] = 0;
    }
}



void relation::read_from_relation(relation* input, int full_delta)
{
    std::vector<u64> prefix = {};
    vector_buffer vb;
    vb.vector_buffer_create_empty();

    full[mcomm.get_rank()]->as_vector_buffer_recursive(&vb, prefix);

    if (full_delta == DELTA)
        populate_delta(vb.size/sizeof(u64), (u64*)vb.buffer);

    else if (full_delta == FULL)
        populate_full(vb.size/sizeof(u64), (u64*)vb.buffer);

    vb.vector_buffer_free();
}
#endif


void relation::initialize_relation(mpi_comm& mcomm)
{
    /// Main : Execute : init : buffer_init : start
    this->mcomm = mcomm;

    u32 buckets = mcomm.get_local_nprocs();

    default_sub_bucket_per_bucket_count = 1;
    int rank = mcomm.get_local_rank();
    int nprocs = mcomm.get_local_nprocs();

    newt_element_count = 0;
    full_element_count = 0;
    delta_element_count = 0;
    delta = new google_relation[buckets];
    full = new google_relation[buckets];
    newt = new google_relation[buckets];

    sub_bucket_per_bucket_count = new u32[buckets];
    for (u32 b = 0; b < buckets; b++)
    {
        //delta[b] = new google_relation();
        //full[b] = new google_relation();
        //newt[b] = new google_relation();

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
    for (u32 b = 0; b < buckets; b++)
    {
        sub_bucket_rank[b] = new u32[sub_bucket_per_bucket_count[b]];
        std::unordered_set<int> distinct_ranks;
        for (u64 x = 0; x < sub_bucket_per_bucket_count[b]; x++)
        {
            sub_bucket_rank[b][x] = rcount % nprocs;
            set_last_rank(rcount);

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

    /// Main : Execute : init : buffer_init : end


    /// reading from file
    if (initailization_type != -1)
    {
        /// Main : Execute : init : io : end

        file_io.parallel_read_input_relation_from_file_to_local_buffer(arity, filename, mcomm.get_local_comm());

        file_io.buffer_data_to_hash_buffer_col(arity, filename, join_column_count, get_bucket_count(), sub_bucket_rank, sub_bucket_per_bucket_count, mcomm.get_local_comm());

        file_io.delete_raw_buffers();

        /* Copy data from buffer to relation */
        if (initailization_type == DELTA)
            populate_delta(file_io.get_hash_buffer_size(), file_io.get_hash_buffer());

        else if (initailization_type == FULL)
            populate_full(file_io.get_hash_buffer_size(), file_io.get_hash_buffer());

        //int f_size = get_full_element_count();
        //u32 g_f_size = 0;
        //MPI_Allreduce(&f_size, &g_f_size, 1, MPI_INT, MPI_SUM, mcomm.get_local_comm());
        //if (rank == 0)
        //    std::cout << "After Init " << filename << " " << g_f_size << std::endl;
        //std::cout << "Writing " << file_io.get_hash_buffer_size() << " bytes" << std::endl;

        file_io.delete_hash_buffers();
    }
}


void relation::initialize_relation_in_scc(bool init_status)
{
    if (init_status == true)
    {
        //if (mcomm.get_local_rank() == 0)
        //    std::cout << "Shift will take place for " << get_debug_id() << std::endl;
        insert_full_in_delta();
    }
}



void relation::populate_full(int buffer_size, u64* buffer)
{
    u32 counter = 0;
    u64 t[arity+1];
    u32 buckets = get_bucket_count();

    for (int i = 0; i < buffer_size; i = i + (arity+1))
    {
        uint64_t bucket_id = tuple_hash(buffer + i, join_column_count) % buckets;

        for (u32 a = i; a < i + arity + 1; a++)
            t[a-i] = buffer[a];

        if (full[bucket_id].insert_tuple_from_array(t, (arity+1)) == true)
        {
            full_element_count++;
            full_bucket_element_count[bucket_id]++;
            counter++;
        }
    }
}



void relation::populate_delta (int buffer_size, u64* buffer)
{
    u64 t[arity+1];
    u32 buckets = get_bucket_count();

    for (int i = 0; i < buffer_size; i = i + (arity+1))
    {
        uint64_t bucket_id = tuple_hash(buffer + i, join_column_count) % buckets;

        for (u32 a = i; a < i + arity + 1; a++)
            t[a-i] = buffer[a];

        if (delta[bucket_id].insert_tuple_from_array(t, arity+1) == true)
        {
            delta_element_count++;
            delta_bucket_element_count[bucket_id]++;
        }
    }
}



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

    for (u32 i = 0; i < buckets; i++)
    {
        full[i].remove_tuple();
        delta[i].remove_tuple();
        newt[i].remove_tuple();

        //delete[] full[i];
        //delete[] delta[i];
        //delete[] newt[i];

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



void relation::copy_relation(relation*& recv_rel, mpi_comm output_comm, int target_cumulative_rank, int tuples_per_task, u32 input_buckets, u32 output_buckets)
{
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

    int full_output_size=0;
    vector_buffer* full_output = new vector_buffer[output_comm.get_nprocs()];
    for (int i = 0; i < output_comm.get_nprocs(); i++)
        full_output[i].vector_buffer_create_empty();

    int delta_output_size=0;
    vector_buffer *delta_input_buffer = new vector_buffer[input_buckets];
    for (u32 i = 0; i < input_buckets; i++)
        full_input_buffer[i].vector_buffer_create_empty();

    int delta_process_size[output_comm.get_nprocs()];
    memset(delta_process_size, 0, output_comm.get_nprocs() * sizeof(int));

    vector_buffer* delta_output = new vector_buffer[output_comm.get_nprocs()];
    for (int i = 0; i < output_comm.get_nprocs(); i++)
        delta_output[i].vector_buffer_create_empty();


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
                for (u32 j =0; j < arity; j++)
                    memcpy(reordered_cur_path + j, (&full_input_buffer[i])->buffer + ((s + j) * sizeof(u64)), sizeof(u64));

                uint64_t bucket_id = tuple_hash(reordered_cur_path, join_column_count) % output_buckets;
                uint64_t sub_bucket_id = tuple_hash(reordered_cur_path + join_column_count, arity-join_column_count) % output_sub_bucket_per_bucket_count[bucket_id];
                int index = output_sub_bucket_rank[bucket_id][sub_bucket_id] + target_cumulative_rank;

                full_process_size[index] = full_process_size[index] + arity;
                full_output[index].vector_buffer_append((const unsigned char*)reordered_cur_path, sizeof(u64)*arity);
                full_output_size = full_output_size + (arity+1);
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
                for (u32 j =0; j < arity; j++)
                    memcpy(reordered_cur_path + j, (&delta_input_buffer[i])->buffer + ((s + j) * sizeof(u64)), sizeof(u64));

                uint64_t bucket_id = tuple_hash(reordered_cur_path, join_column_count) % output_buckets;
                uint64_t sub_bucket_id = tuple_hash(reordered_cur_path + join_column_count, arity-join_column_count) % output_sub_bucket_per_bucket_count[bucket_id];
                int index = output_sub_bucket_rank[bucket_id][sub_bucket_id] + target_cumulative_rank;

                delta_process_size[index] = delta_process_size[index] + arity;
                delta_output[index].vector_buffer_append((const unsigned char*)reordered_cur_path, sizeof(u64)*arity);
                delta_output_size = delta_output_size + (arity+1);
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


    int full_buffer_size;
    u64* full_buffer;
    all_to_all_comm(full_output, full_output_size, full_process_size, &full_buffer_size, &full_buffer, mcomm.get_comm());
    delete[] full_output;


    finalize_relation();
    recv_rel->set_initailization_type(-1);
    recv_rel->initialize_relation(output_comm);


    u64 t[arity];
    for (int i=0; i < full_buffer_size; i=i+arity)
    {
        for (u32 j = 0; j < arity; j++)
            t[j] = full_buffer[i+j];
        recv_rel->insert_in_full(t);
    }


    int delta_buffer_size;
    u64* delta_buffer;
    all_to_all_comm(delta_output, delta_output_size, delta_process_size, &delta_buffer_size, &delta_buffer, mcomm.get_comm());
    delete[] delta_output;



    for (int i=0; i < delta_buffer_size; i=i+arity)
    {
        for (u32 j = 0; j < arity; j++)
            t[j] = full_buffer[i+j];
        recv_rel->insert_in_delta(t);
    }

    delete[] full_buffer;
    delete[] delta_buffer;
}



bool relation::find_in_full(u64* t, int length)
{
    uint64_t bucket_id = tuple_hash(t, join_column_count) % get_bucket_count();
    return full[bucket_id].find_tuple_from_array(t, length);
}



bool relation::find_in_delta(u64* t, int length)
{
    uint64_t bucket_id = tuple_hash(t, join_column_count) % get_bucket_count();
    return delta[bucket_id].find_tuple_from_array(t, length);
}



bool relation::find_in_newt(u64* t, int length)
{
    uint64_t bucket_id = tuple_hash(t, join_column_count) % get_bucket_count();
    return newt[bucket_id].find_tuple_from_array(t, length);
}



bool relation::insert_in_delta(u64* t)
{
    uint64_t bucket_id = tuple_hash(t, join_column_count) % get_bucket_count();
    u32 sub_bucket_id = 0;
    if (is_canonical == false)
        sub_bucket_id = tuple_hash(t + join_column_count, arity-join_column_count) % sub_bucket_per_bucket_count[bucket_id];

    assert((int)bucket_id == mcomm.get_local_rank());
    if (delta[bucket_id].insert_tuple_from_array(t, arity+1) == true)
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
    uint64_t bucket_id = tuple_hash(t, join_column_count) % get_bucket_count();
    u32 sub_bucket_id = 0;
    if (is_canonical == false)
        sub_bucket_id = tuple_hash(t + join_column_count, arity-join_column_count) % sub_bucket_per_bucket_count[bucket_id];

    if((int)bucket_id != mcomm.get_local_rank())
    {
        for (u32 x=0; x < arity+1; x++)
            std::cout << t[x] << " ";
        std::cout << "BID " << bucket_id << " sub-bucket id " << sub_bucket_id << " rank " << mcomm.get_local_rank() << std::endl;
    }
    if (newt[bucket_id].insert_tuple_from_array(t, arity+1) == true)
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
    uint64_t bucket_id = tuple_hash(t, join_column_count) % buckets;
    uint64_t sub_bucket_id=0;
    if (is_canonical == false)
        sub_bucket_id = tuple_hash(t + join_column_count, arity-join_column_count) % sub_bucket_per_bucket_count[bucket_id];

    assert((int)bucket_id == mcomm.get_local_rank());

#if 0
    if (get_debug_id() == "rel_path_2_1")
    {
    std::cout << "Tuples to insert " << full << " ";
    for (u32 i=0; i < arity+1; i++)
        std::cout << t[i] << " ";
    std::cout << std::endl;
    }
#endif

    if (full[bucket_id].insert_tuple_from_array(t, arity+1) == true)
    {
        full_element_count++;
        full_bucket_element_count[bucket_id]++;
        full_sub_bucket_element_count[bucket_id][sub_bucket_id]++;
        bucket_map[bucket_id] = 1;

        //std::cout << "Full size " << full_element_count << std::endl;
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
    {
        input_buffer[i].vector_buffer_create_empty();
        if (bucket_map[i] == 1)
        {
            std::vector<u64> prefix = {};
            delta[i].as_vector_buffer_recursive(&(input_buffer[i]), prefix);
            for (u64 j = 0; j < (&input_buffer[i])->size / sizeof(u64); j=j+(arity+1))
            {
                if (insert_in_full ( (u64*)( (input_buffer[i].buffer) + (j*sizeof(u64)) )) == true)
                    insert_success++;
            }
            delta[i].remove_tuple();
            input_buffer[i].vector_buffer_free();
        }
    }
    //if (get_debug_id() == "rel_edge_2_2" || get_debug_id() == "rel_edge_2_1_2")
    //if (mcomm.get_local_rank() == 0)
    //    std::cout << "[" << get_debug_id() << "] inserting delta in full insert_success " << insert_success << std::endl;

    set_delta_element_count(0);
    delete[] input_buffer;

    return insert_success;
}



int relation::insert_full_in_delta()
{
    u32 insert_success = 0;
    u32 buckets = get_bucket_count();
    vector_buffer *input_buffer = new vector_buffer[buckets];

    for (u32 i = 0; i < buckets; i++)
    {
        input_buffer[i].vector_buffer_create_empty();
        if (bucket_map[i] == 1)
        {
            std::vector<u64> prefix = {};
            full[i].as_vector_buffer_recursive(&(input_buffer[i]), prefix);
            //std::cout << debug_id << " SIZE OF BUFFER " << (&input_buffer[i])->size << std::endl;
            for (u64 j = 0; j < (&input_buffer[i])->size / sizeof(u64); j=j+(arity+1))
            {
                if (insert_in_delta ( (u64*)( (input_buffer[i].buffer) + (j*sizeof(u64)) )) == true)
                    insert_success++;
            }
            full[i].remove_tuple();
            input_buffer[i].vector_buffer_free();
        }
    }

    set_full_element_count(0);
    delete[] input_buffer;

    return insert_success;
}



void relation::local_insert_in_delta()
{
    int rank;
    MPI_Comm_rank(mcomm.get_comm(), &rank);
    u32 buckets = get_bucket_count();

    delete[] delta;


    delta = newt;

    /*
    u32 i = mcomm.get_rank();
    vector_buffer *vb_delta = new vector_buffer[buckets];
    vb_delta[i].vector_buffer_create_empty();
    std::vector<u64> prefix = {};
    delta[i].as_vector_buffer_recursive(&(vb_delta[i]), prefix);

    //if (i == 0)
    //    std::cout << "[" << get_debug_id() << "] Test " << mcomm.get_rank() << " DELTA " << vb_delta[i].size/(sizeof(u64) * (arity + 1)) << " arity " << arity + 1 << std::endl;

    vb_delta[i].vector_buffer_free();

    delete[] vb_delta;
    */


    delta_element_count = newt_element_count;
    //if (rank == 0)
    //    std::cout << "[" << get_debug_id() << "] copyng newt pointer to delta   " << delta_element_count << std::endl;

    memcpy(delta_bucket_element_count, newt_bucket_element_count, buckets * sizeof(u32));
    for (u32 b = 0; b < buckets; b++)
    {
        memcpy(delta_sub_bucket_element_count[b], newt_sub_bucket_element_count[b], sub_bucket_per_bucket_count[b] * sizeof(u32));
        memset(newt_sub_bucket_element_count[b], 0, sub_bucket_per_bucket_count[b] * sizeof(u32));
    }

    newt = new google_relation[buckets];
    //for(u32 i=0; i<buckets; i++)
    //    newt[i] = new google_relation();
    newt_element_count = 0;
    memset(newt_bucket_element_count, 0, buckets * sizeof(u32));
}
