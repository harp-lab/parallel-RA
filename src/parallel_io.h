#ifndef PARALLEL_H
#define PARALLEL_H

#include <mpi.h>



inline u64 tunedhash(const u8* bp, const u32 len)
{
    u64 h0 = 0xb97a19cb491c291d;
    u64 h1 = 0xc18292e6c9371a17;
    const u8* const ep = bp+len;
    while (bp < ep)
    {
        h1 ^= *bp;
        h1 *= 31;
        h0 ^= (((u64)*bp) << 17) ^ *bp;
        h0 *= 0x100000001b3;
        h0 = (h0 >> 7) | (h0 << 57);
        ++bp;
    }

    return h0 ^ h1;// ^ (h1 << 31);
}



u64 outer_hash(const u64 val)
{
    return tunedhash((u8*)(&val),sizeof(u64));
}


class parallel_io
{

private:

    // filename of the data
    const char *file_name;

    // defaulted to 2
    u32 col_count;

    // total number of rows / nprocs
    u32 entry_count;
    u64 *input_buffer;

    // total number of rows after hashing
    u64* hash_buffer;
    u64 hash_buffer_size;

public:


    parallel_io(const char *fname)
    {
        file_name = fname;
    }


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
        return hash_buffer_size;
    }

    void parallel_read_input_relation_from_file_to_local_buffer(int rank, int nprocs)
    {
        u32 global_row_count;
        u32 global_col_count;

        if (rank == 0)
        {
            char meta_data_filename[1024];
            sprintf(meta_data_filename, "%s/meta_data.txt", file_name);
            printf("Opening File %s\n", meta_data_filename);

            FILE *fp_in;
            fp_in = fopen(meta_data_filename, "r");
            if (fscanf (fp_in, "(row count)\n%d\n(col count)\n%d", &global_row_count, &global_col_count) != 2)
            {
                printf("Wrong input format (Meta Data)\n");
                MPI_Abort(MPI_COMM_WORLD, -1);
            }
            fclose(fp_in);
        }

        MPI_Bcast(&global_row_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&global_col_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

        col_count = global_col_count;

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

        input_buffer = new u64[entry_count * global_col_count];
        u32 rb_size = pread(fp, input_buffer, entry_count * global_col_count * sizeof(u64), read_offset * global_col_count * sizeof(u64));
        if (rb_size != entry_count * global_col_count * sizeof(u64))
        {
            std::cout << "Wrong IO: rank: " << rank << " " << rb_size << " " <<  entry_count << " " << global_col_count << std::endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
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

    void buffer_data_to_hash_buffer(int hash_column_index, u32 nprocs, MPI_Comm comm)
    {
        /* process_size[j] stores the number of samples to be sent to process with rank j */
        int* process_size = new int[nprocs];
        memset(process_size, 0, nprocs * sizeof(int));

        /* vector[i] contains the data that needs to be sent to process i */
        std::vector<u64> process_data_vector[nprocs];

        printf("entry_count = %d nprocs = %d\n", (int)entry_count, (int)nprocs);

        if (hash_column_index == 0)
        {
            for (u32 i = 0; i < entry_count * col_count; i=i+2)
            {
                uint64_t index = outer_hash(input_buffer[i])%nprocs;
                process_size[index] = process_size[index] + col_count;

                process_data_vector[index].push_back(input_buffer[i]);
                process_data_vector[index].push_back(input_buffer[i + 1]);

                //if (rank == 0)
                    //std::cout << input_buffer[i] << " : " << input_buffer[i + 1] << std::endl;
            }
        }
        else
        {
            for (u32 i = 0; i < entry_count * col_count; i=i+2)
            {
                uint64_t index = outer_hash(input_buffer[i+1])%nprocs;
                process_size[index] = process_size[index] + col_count;

                process_data_vector[index].push_back(input_buffer[i]);
                process_data_vector[index].push_back(input_buffer[i + 1]);
            }
        }


        int prefix_sum_process_size[nprocs];
        memset(prefix_sum_process_size, 0, nprocs * sizeof(int));
        for (u32 i = 1; i < nprocs; i++)
            prefix_sum_process_size[i] = prefix_sum_process_size[i - 1] + process_size[i - 1];


        int process_data_buffer_size = prefix_sum_process_size[nprocs - 1] + process_size[nprocs - 1];


        u64* process_data = new u64[process_data_buffer_size];

        for (u32 i = 0; i < nprocs; i++)
            memcpy(process_data + prefix_sum_process_size[i], &process_data_vector[i][0], process_data_vector[i].size() * sizeof(u64));


        /* This step prepares for actual data transfer */
        /* Every process sends to every other process the amount of data it is going to send */
        int recv_process_size_buffer[nprocs];
        memset(recv_process_size_buffer, 0, nprocs * sizeof(int));
        MPI_Alltoall(process_size, 1, MPI_INT, recv_process_size_buffer, 1, MPI_INT, comm);


        int prefix_sum_recv_process_size_buffer[nprocs];
        memset(prefix_sum_recv_process_size_buffer, 0, nprocs * sizeof(int));
        for (u32 i = 1; i < nprocs; i++)
            prefix_sum_recv_process_size_buffer[i] = prefix_sum_recv_process_size_buffer[i - 1] + recv_process_size_buffer[i - 1];


        /* Sending data to all processes */
        /* What is the buffer size to allocate */
        hash_buffer_size = 0;
        for(u32 i = 0; i < nprocs; i++)
            hash_buffer_size = hash_buffer_size + recv_process_size_buffer[i];


        printf("hash_buffer_size = %d\n", (int)hash_buffer_size);
        hash_buffer = new u64[hash_buffer_size];
        MPI_Alltoallv(process_data, process_size, prefix_sum_process_size, MPI_UNSIGNED_LONG_LONG, hash_buffer, recv_process_size_buffer, prefix_sum_recv_process_size_buffer, MPI_UNSIGNED_LONG_LONG, comm);

        hash_buffer_size = hash_buffer_size / col_count;

        delete[] process_data;
        delete[] process_size;

        return;
    }
};

#endif
