/*
 * Parallel IO of relations
 * Copyright (c) Sidharth Kumar, et al, see License.md
 */



#ifndef _PARALLEL_IO_H
#define _PARALLEL_IO_H


class parallel_io
{

private:

    /// filename of the data
    const char *file_name;

    /// arity of the relation
    int col_count;

    /// total number of rows / nprocs
    int entry_count;
    u64* input_buffer;

    /// total number of rows after hashing
    int hash_buffer_size;
    u64* hash_buffer;

public:

    parallel_io();

    u64* get_hash_buffer()  {  return hash_buffer;  }
    int get_hash_buffer_size()  {  return hash_buffer_size;  }


    void delete_raw_buffers()   {   if (entry_count != 0)        delete[] input_buffer;}
    void delete_hash_buffers()  {  delete[] hash_buffer;  }


    /// read file with offset
    void parallel_read_input_relation_from_file_with_offset(u32 arity, const char *fname, MPI_Comm lcomm);


    /// offset reads (parallel IO)
    void parallel_read_input_relation_from_file_to_local_buffer(u32 arity, const char *fname, MPI_Comm lcomm);


    /// move tuples to the appropriate process
    void buffer_data_to_hash_buffer_col(u32 arity, const char *fname, u32 join_column_count, u32 buckets, u32** sub_bucket_rank, u32* sub_bucket_count, MPI_Comm lcomm);
};

#endif
