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

#include "btree.h"
#include "btree_relation.h"

static int rank = 0;
static u32 nprocs = 1;

void all_to_all_test(u64 data_size, u32 iteration_index)
{
    double c1 = 0;
    double i2 = 0;

    int send_process_size[nprocs];
    int recv_process_size[nprocs];

    int send_process_prefix[nprocs];
    int recv_process_prefix[nprocs];
    send_process_prefix[0] = 0;
    recv_process_prefix[0] = 0;

    send_process_size[0] = data_size / nprocs;
    recv_process_size[0] = data_size / nprocs;

    for(u32 i = 1; i < nprocs; i++)
    {
        send_process_prefix[i] = send_process_prefix[i - 1] +  send_process_size[i - 1];
        send_process_size[i] = data_size / nprocs;

        recv_process_prefix[i] = recv_process_prefix[i - 1] + recv_process_size[i - 1];
        recv_process_size[i] = data_size / nprocs;
    }

    u64* send_buffer = 0;
    send_buffer = new u64[data_size];
    memset(send_buffer, 0, data_size * sizeof(u64));

    for(u32 i = 0; i < data_size; i++)
        send_buffer[i] = i / (data_size/nprocs);


    u64 *recv_buffer = 0;
    recv_buffer = new u64[data_size];
    memset(recv_buffer, 0, data_size * sizeof(u64));

    c1 = MPI_Wtime();
    MPI_Alltoallv(send_buffer, send_process_size, send_process_prefix, MPI_UNSIGNED_LONG_LONG, recv_buffer, recv_process_size, recv_process_prefix, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);

    i2 = MPI_Wtime();

    if (rank == 0)
    for (u32 i = 0; i < data_size; i++)
        if (recv_buffer[i] != (u64)rank)
        {
            std::cout << "Error !!!!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

    delete[] send_buffer;
    delete[] recv_buffer;

    double time = i2 - c1;
    double max_time;
    MPI_Allreduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (time == max_time)
    std::cout << "[" << iteration_index << "] " << nprocs << " " << data_size << " "  << " Rank " << rank <<  " All to all communication time: " << max_time << " Rank 0 time: " << (i2 - c1) << std::endl;


}


int main(int argc, char **argv)
{
    // Initializing MPI
    MPI_Init(&argc, &argv);
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    nprocs = size;

    for (int i = 33; i <= 36; i++)
    {
        all_to_all_test(pow(2, i)/nprocs, i);
        MPI_Barrier(MPI_COMM_WORLD);
    }


    for (int i = 33; i <= 36; i++)
    {
        all_to_all_test(pow(2, i)/nprocs, i);
        MPI_Barrier(MPI_COMM_WORLD);
    }


    // Finalizing MPI
    MPI_Finalize();

    return 0;
}








