//#include <chrono>
#include <sys/stat.h>
#include <errno.h>
#include <limits.h>
#include <arpa/inet.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <iostream>
#include <fstream>


#include "btree.h"
#include "btree_relation.h"

#define COL_COUNT 2

static int rank = 0;
static u32 nprocs = 1;
static MPI_Comm comm;


#if 1
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

    return h0 ^ h1;
}



u64 outer_hash(const u64 val)
{
    return tunedhash((u8*)(&val),sizeof(u64));
}
#endif


void parallel_read_input_relation_from_file_to_local_buffer(const char *file_name, u64** read_buffer, u64 row_count, u64* local_row_count)
{
    u64 read_offset;
    read_offset = ceil((float)row_count / nprocs) * rank;

    if (read_offset > row_count)
    {
        *local_row_count = 0;
        return;
    }

    if (read_offset + ceil((float)row_count / nprocs) > row_count)
      *local_row_count = row_count - read_offset;
    else
      *local_row_count = (u64) ceil((float)row_count / nprocs);

    if (*local_row_count == 0)
        return;

    char data_filename[1024];
    sprintf(data_filename, "%s/data.raw", file_name);
    int fp = open(data_filename, O_RDONLY);


    size_t bsize = (u64) *local_row_count *  (u64)COL_COUNT * sizeof(u64);
    *read_buffer = new u64[*local_row_count * COL_COUNT];


    u64 rb_size = pread(fp, *read_buffer, bsize, read_offset * COL_COUNT * sizeof(u64));
    if (rb_size != *local_row_count * COL_COUNT * sizeof(u64))
    {
        std::cout << "Wrong IO: rank: " << rank << " " << rb_size << " " <<  bsize << " " << *local_row_count << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    close(fp);


    /*
    std::ifstream myFile (data_filename, std::ios::in | std::ios::binary);
    myFile.seekg (read_offset * COL_COUNT * sizeof(u64));
    myFile.read ((char*)*read_buffer, bsize);
    myFile.close();
    */


    return;
}

void buffer_data_to_hash_buffer(u64 local_number_of_rows, u64* input_data, u64** outer_hash_data, u64* outer_hash_buffer_size, MPI_Comm comm, u32 index)
{
    /* process_size[j] stores the number of samples to be sent to process with rank j */
    int process_size[nprocs];
    memset(process_size, 0, nprocs * sizeof(int));

    /* vector[i] contains the data that needs to be sent to process i */
    std::vector<tuple<2>> *process_data_vector;
    process_data_vector = new std::vector<tuple<2>>[nprocs];

    tuple<2> dt;
    for (u32 i = 0; i < local_number_of_rows * COL_COUNT; i=i+2)
    {
        uint64_t index = outer_hash(input_data[i])%nprocs;
        process_size[index] = process_size[index] + COL_COUNT;

        dt[0] = input_data[i];
        dt[1] = input_data[i + 1];

        process_data_vector[index].push_back(dt);
    }

    int prefix_sum_process_size[nprocs];
    memset(prefix_sum_process_size, 0, nprocs * sizeof(int));

    process_size[0] = 2 * process_data_vector[0].size();
    if (rank == 0 && index == 1)
        std::cout << "PS 0" << " " << process_size[0] << std::endl;

    for (u32 i = 1; i < nprocs; i++)
    {
        prefix_sum_process_size[i] = prefix_sum_process_size[i - 1] + process_size[i - 1];
        process_size[i] = 2 * process_data_vector[i].size();

        if (rank == 0 && index == 1)
            std::cout << "PS " << i << " " << process_size[i] << std::endl;
    }


    int process_data_buffer_size = prefix_sum_process_size[nprocs - 1] + process_size[nprocs - 1];
    u64* process_data = new u64[process_data_buffer_size];

    for (u32 i = 0; i < nprocs; i++)
        memcpy(process_data + prefix_sum_process_size[i], &process_data_vector[i][0], process_data_vector[i].size() * sizeof(tuple<2>));

    delete[] process_data_vector;

    /* This step prepares for actual data transfer */
    /* Every process sends to every other process the amount of data it is going to send */
    int recv_process_size_buffer[nprocs];
    memset(recv_process_size_buffer, 0, nprocs * sizeof(int));
    MPI_Alltoall(process_size, 1, MPI_INT, recv_process_size_buffer, 1, MPI_INT, comm);


    /* Sending data to all processes, what is the buffer size to allocate */
    int prefix_sum_recv_process_size_buffer[nprocs];
    memset(prefix_sum_recv_process_size_buffer, 0, nprocs * sizeof(int));
    *outer_hash_buffer_size = recv_process_size_buffer[0];
    if (rank == 0 && index == 1)
        std::cout << "OHBS 0" << " " << *outer_hash_buffer_size << std::endl;
    for (u32 i = 1; i < nprocs; i++)
    {
        prefix_sum_recv_process_size_buffer[i] = prefix_sum_recv_process_size_buffer[i - 1] + recv_process_size_buffer[i - 1];
        *outer_hash_buffer_size = *outer_hash_buffer_size + recv_process_size_buffer[i];

        if (rank == 0 && index == 1)
            std::cout << "OHBS " << i << " " << *outer_hash_buffer_size << std::endl;
    }

    *outer_hash_data = new u64[*outer_hash_buffer_size];
    MPI_Alltoallv(process_data, process_size, prefix_sum_process_size, MPI_UNSIGNED_LONG_LONG, *outer_hash_data, recv_process_size_buffer, prefix_sum_recv_process_size_buffer, MPI_UNSIGNED_LONG_LONG, comm);

    *outer_hash_buffer_size = *outer_hash_buffer_size / COL_COUNT;

    delete[] process_data;

    return;
}


int main(int argc, char **argv)
{
    // Initializing MPI
    MPI_Init(&argc, &argv);
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    comm = MPI_COMM_WORLD;

    nprocs = (u32)size;

    u32 relation_count = atoi(argv[1]);
    u64 *  row_count = new u64[relation_count];
    u64 *  hash_entry_count = new u64[relation_count];

    u64 *  local_entry_count = new u64[relation_count];
    u64 ** input_buffer = new u64*[relation_count];
    u64 ** hashed_data = new u64*[relation_count];

    double * ior_start = new double[relation_count];
    double * ior_end = new double[relation_count];
    double * hash_start = new double[relation_count];
    double * hash_end = new double[relation_count];
    double * union_start = new double[relation_count];
    double * union_end = new double[relation_count];


    relation<2> unionG;
    tuple<2> t1;

    double start = MPI_Wtime();
    u64 total_tuple = 0;
    u64 inserted_tuple = 0;
    for (u32 i = 0; i < relation_count; i++)
    {
        ior_start[i] = MPI_Wtime();
        row_count[i] = atoi(argv[3 + 2*i]);
        parallel_read_input_relation_from_file_to_local_buffer(argv[2 + 2*i], &(input_buffer[i]), row_count[i], &local_entry_count[i]);
        ior_end[i] = MPI_Wtime();

        hash_start[i] = MPI_Wtime();
        buffer_data_to_hash_buffer(local_entry_count[i], input_buffer[i], &(hashed_data[i]), &(hash_entry_count[i]), MPI_COMM_WORLD, i);
        hash_end[i] = MPI_Wtime();

        union_start[i] = MPI_Wtime();
        for ( u32 j = 0; j < hash_entry_count[i] * COL_COUNT; j = j + 2)
        {
            total_tuple++;
            t1[0] = hashed_data[i][j];
            t1[1] = hashed_data[i][j+1];
            if (unionG.insert(t1))
                inserted_tuple++;
        }
        union_end[i] = MPI_Wtime();
    }
    double end = MPI_Wtime();

    double elapsed_time = end - start;
    double max_time = 0;
    u64 total_sum = 0;
    u64 total_sum2 = 0;
    MPI_Allreduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&inserted_tuple, &total_sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
    MPI_Allreduce(&total_tuple, &total_sum2, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);

    if (max_time == elapsed_time)
    {
        std::cout << nprocs << " : "
                  << "[L] Total tuple: " << total_tuple
                  << " [G] Total tuple: " << total_sum2
                  << " [L] Inserted tuple: " << inserted_tuple
                  << " [G] Inserted tuple: " << total_sum
                  << std::endl;

        double time = 0;
        for (u32 i = 0; i < relation_count; i++)
        {
            std::cout << "[" << i << "]: "
                      << " Global row count: " << row_count[i]
                      << " Local row count: " << local_entry_count[i]
                      << " Hash row count: " << hash_entry_count[i]
                      << " Read time: " << (ior_end[i] - ior_start[i])
                      << " Hash time: " << (hash_end[i] - hash_start[i])
                      << " Insert time: " << (union_end[i] - union_start[i])
                      << std::endl;

            time = time + (ior_end[i] - ior_start[i]) + (hash_end[i] - hash_start[i]) + (union_end[i] - union_start[i]);
        }

        std::cout << "Read + hash + insert: " << time << std::endl;
        std::cout << "Total time: " << max_time << " " << time << std::endl;
    }

    /*

    for (u32 i = 0; i < relation_count; i++)
    {
        delete[] input_buffer[i];
        delete[] hashed_data[i];
    }

    delete[] input_buffer;
    delete[] hashed_data;


    delete[]  row_count;
    delete[]  hash_entry_count;
    delete[]  local_entry_count;


    delete[] ior_start;
    delete[] ior_end;
    delete[] hash_start;
    delete[] hash_end;
    delete[] union_start;
    delete[] union_end;
    */

    // Finalizing MPI
    MPI_Finalize();

    return 0;
}


// 256 : [L] Total tuple: 99583402 [G] Total tuple: 664659334 [L] Inserted tuple: 99583240 [G] Inserted tuple: 664616755 Write time: 24.5067
//4096 : [L] Total tuple: 97510090 [G] Total tuple: 664659334 [L] Inserted tuple: 97510078 [G] Inserted tuple: 664616755 Write time: 21.0971







