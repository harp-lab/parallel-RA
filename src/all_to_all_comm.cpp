#include "parallel_RA_inc.h"


void all_to_all_comm(vector_buffer* local_join_output, int* process_size, u64 *outer_hash_buffer_size, u64 **outer_hash_data, MPI_Comm comm)
{
    int nprocs;
    MPI_Comm_size(comm, &nprocs);

    /* This step prepares for actual data transfer */
    /* Every process sends to every other process the amount of data it is going to send */
    int recv_process_size_buffer[nprocs];
    memset(recv_process_size_buffer, 0, nprocs * sizeof(int));
    MPI_Alltoall(process_size, 1, MPI_INT, recv_process_size_buffer, 1, MPI_INT, comm);

    int prefix_sum_process_size[nprocs];
    memset(prefix_sum_process_size, 0, nprocs * sizeof(int));

    for (int i = 1; i < nprocs; i++)
        prefix_sum_process_size[i] = prefix_sum_process_size[i - 1] + process_size[i - 1];

    int process_data_buffer_size = prefix_sum_process_size[nprocs - 1] + process_size[nprocs - 1];

    u64* process_data = 0;
    process_data = new u64[process_data_buffer_size];
    memset(process_data, 0, process_data_buffer_size * sizeof(u64));

    for(int i = 0; i < nprocs; i++)
    {
        memcpy(process_data + prefix_sum_process_size[i], (&local_join_output[i])->buffer, (&local_join_output[i])->size);
        local_join_output[i].vector_buffer_free();
    }

    int prefix_sum_recv_process_size_buffer[nprocs];
    memset(prefix_sum_recv_process_size_buffer, 0, nprocs * sizeof(int));

    /* Sending data to all processes: What is the buffer size to allocate */
    *outer_hash_buffer_size = recv_process_size_buffer[0];
    for(int i = 1; i < nprocs; i++)
    {
        prefix_sum_recv_process_size_buffer[i] = prefix_sum_recv_process_size_buffer[i - 1] + recv_process_size_buffer[i - 1];
        *outer_hash_buffer_size = *outer_hash_buffer_size + recv_process_size_buffer[i];
    }

    *outer_hash_data = new u64[*outer_hash_buffer_size];
    memset(*outer_hash_data, 0, *outer_hash_buffer_size * sizeof(u64));


    MPI_Alltoallv(process_data, process_size, prefix_sum_process_size, MPI_UNSIGNED_LONG_LONG, *outer_hash_data, recv_process_size_buffer, prefix_sum_recv_process_size_buffer, MPI_UNSIGNED_LONG_LONG, comm);

    delete[] process_data;

    return;
}
