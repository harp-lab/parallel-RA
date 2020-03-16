#ifndef BALANCED_PARALLEL_H
#define BALANCED_PARALLEL_H

#include "vector_buffer.h"
#include <mpi.h>



class parallel_io
{

private:

    int rank;
    int nprocs;
    MPI_Comm world_comm;

    /// filename of the data
    const char *file_name;

    /// defaulted to 2
    u32 col_count;

    /// total number of rows / nprocs
    u32 entry_count;
    u64* input_buffer;

    /// total number of rows after hashing
    u64 hash_buffer_size;
    u64* hash_buffer;

public:


    u64* get_hash_buffer()
    {
        return hash_buffer;
    }

    u64 get_col_count()
    {
        return col_count;
    }

    u64 get_row_count()
    {
        return entry_count;
    }

    u64 get_hash_buffer_size()
    {
        return hash_buffer_size;
    }

    void parallel_read_input_relation_from_file_to_local_buffer(const char *fname, MPI_Comm comm)
    {
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &nprocs);
        world_comm = comm;
        file_name = fname;

        u32 global_row_count;

        /* Read the metadata file containing the total number of rows and columns */
        if (rank == 0)
        {
            char meta_data_filename[1024];
            sprintf(meta_data_filename, "%s/meta_data.txt", file_name);
            printf("Opening File %s\n", meta_data_filename);

            FILE *fp_in;
            fp_in = fopen(meta_data_filename, "r");
            if (fscanf (fp_in, "(row count)\n%d\n(col count)\n%d", &global_row_count, &col_count) != 2)
            {
                printf("Wrong input format (Meta Data)\n");
                MPI_Abort(comm, -1);
            }
            fclose(fp_in);
        }

        /* Broadcast the total number of rows and column to all processes */
        MPI_Bcast(&global_row_count, 1, MPI_INT, 0, comm);
        MPI_Bcast(&col_count, 1, MPI_INT, 0, comm);


        /* Read all data in parallel */
        u32 read_offset;
        read_offset = ceil((float)global_row_count / nprocs) * rank;

        if (read_offset > global_row_count)
        {
            entry_count = 0;
            return;
        }

        if (read_offset + ceil((float)global_row_count / nprocs) > global_row_count)
            entry_count = global_row_count - read_offset;
        else
            entry_count = (u32) ceil((float)global_row_count / nprocs);

        if (entry_count == 0)
            return;

        char data_filename[1024];
        sprintf(data_filename, "%s/data.raw", file_name);
        int fp = open(data_filename, O_RDONLY);

        input_buffer = new u64[entry_count * col_count];
        u32 rb_size = pread(fp, input_buffer, entry_count * col_count * sizeof(u64), read_offset * col_count * sizeof(u64));
        if (rb_size != entry_count * col_count * sizeof(u64))
        {
            std::cout << "Wrong IO: rank: " << rank << " " << rb_size << " " <<  entry_count << " " << col_count << std::endl;
            MPI_Abort(comm, -1);
        }
        close(fp);

        return;
    }



    void delete_raw_buffers()
    {
        delete[] input_buffer;
    }



    void delete_hash_buffers()
    {
        delete[] hash_buffer;
    }



    /// Hashing based on the first and the second column
    /// The first column decises the bucket id
    /// The second column decides the sub-bucket id
    void buffer_data_to_hash_buffer_col(u32 buckets, u32** sub_bucket_rank, u32* sub_bucket_count)
    {
        /* process_size[i] stores the number of samples to be sent to process with rank i */
        int* process_size = new int[nprocs];
        memset(process_size, 0, nprocs * sizeof(int));

        /* process_data_vector[i] contains the data that needs to be sent to process i */
        vector_buffer* process_data_vector = (vector_buffer*)malloc(sizeof(vector_buffer) * nprocs);
        for (u32 i = 0; i < (u32)nprocs; ++i)
            process_data_vector[i] = vector_buffer_create_empty();

        /* Hashing and buffering data for all to all comm */
        u64 val[col_count];
        for (u32 i = 0; i < entry_count * col_count; i=i+col_count)
        {
            uint64_t bucket_id = hash_function(input_buffer[i]) % buckets;
            uint64_t sub_bucket_id = hash_function(input_buffer[i+1]) % sub_bucket_count[bucket_id];

            int index = sub_bucket_rank[bucket_id][sub_bucket_id];
            process_size[index] = process_size[index] + col_count;

            for (u32 j = 0; j < col_count; j++)
                val[j] = input_buffer[i + j];

            vector_buffer_append(&process_data_vector[index],(unsigned char *) val, sizeof(u64)*col_count);
        }


        /* Transmit the packaged data process_data_vector to all processes */
        all_to_all_comm(process_data_vector, process_size);

        /* Free the data buffer after all to all */
        free (process_data_vector);

        /* Free the process size buffer */
        delete[] process_size;

        return;
    }



    /// All to all communication
    /// process_data_vector has packaged all the data that needs to be sent out
    /// process_size collects the size of the buffer that needs to be sent
    void all_to_all_comm(vector_buffer* process_data_vector, int* process_size)
    {
        /* prefix sum on the send side (required for all to all communication) */
        int prefix_sum_process_size[nprocs];
        memset(prefix_sum_process_size, 0, nprocs * sizeof(int));
        for (u32 i = 1; i < (u32)nprocs; i++)
            prefix_sum_process_size[i] = prefix_sum_process_size[i - 1] + process_size[i - 1];

        /* Total data sent out during all to all */
        int process_data_buffer_size = prefix_sum_process_size[nprocs - 1] + process_size[nprocs - 1];

        /* Buffer to be transmitted during all to all */
        u64* process_data = new u64[process_data_buffer_size];

        /* Populating the buffer to be transmitted */
        for (u32 i = 0; i < (u32)nprocs; i++)
        {
            memcpy(process_data + prefix_sum_process_size[i], (&process_data_vector[i])->buffer, (&process_data_vector[i])->size);
            vector_buffer_free(&process_data_vector[i]);
        }

        /* This step prepares for actual data transfer */
        /* Every process sends to every other process the amount of data it is going to send */
        int recv_process_size_buffer[nprocs];
        memset(recv_process_size_buffer, 0, nprocs * sizeof(int));
        MPI_Alltoall(process_size, 1, MPI_INT, recv_process_size_buffer, 1, MPI_INT, world_comm);


        /* Prefix sum on the receive side (required for all to all communication) */
        int prefix_sum_recv_process_size_buffer[nprocs];
        memset(prefix_sum_recv_process_size_buffer, 0, nprocs * sizeof(int));
        for (u32 i = 1; i < (u32)nprocs; i++)
            prefix_sum_recv_process_size_buffer[i] = prefix_sum_recv_process_size_buffer[i - 1] + recv_process_size_buffer[i - 1];

        /* Total data received during all to all */
        hash_buffer_size = prefix_sum_recv_process_size_buffer[nprocs - 1] + recv_process_size_buffer[nprocs - 1];

        /* Buffer to be reveived after all to all */
        hash_buffer = new u64[hash_buffer_size];

        /* All to all communication */
        MPI_Alltoallv(process_data, process_size, prefix_sum_process_size, MPI_UNSIGNED_LONG_LONG, hash_buffer, recv_process_size_buffer, prefix_sum_recv_process_size_buffer, MPI_UNSIGNED_LONG_LONG, world_comm);

        /* Delete the send buffer */
        delete[] process_data;
    }

};

#endif
