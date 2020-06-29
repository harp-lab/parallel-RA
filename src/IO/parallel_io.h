#ifndef _PARALLEL_IO_H
#define _PARALLEL_IO_H


class parallel_io
{

private:

    /// filename of the data
    const char *file_name;

    /// arity of the relation
    u32 col_count;

    /// total number of rows / nprocs
    u32 entry_count;
    u64* input_buffer;

    /// total number of rows after hashing
    u64 hash_buffer_size;
    u64* hash_buffer;

public:

    u64* get_hash_buffer()  {  return hash_buffer;  }
    u64 get_col_count()  {  return col_count;  }
    u64 get_row_count()  {  return entry_count;  }
    u64 get_hash_buffer_size()  {  return hash_buffer_size;  }
    void parallel_read_input_relation_from_file_to_local_buffer(const char *fname, MPI_Comm lcomm);
    void delete_raw_buffers();
    void delete_hash_buffers()  {  delete[] hash_buffer;  }

    void buffer_data_to_hash_buffer_col(u32 arity, u32 join_column_count, u32 buckets, u32** sub_bucket_rank, u32* sub_bucket_count, MPI_Comm lcomm);
};

#endif
