#include "../src/parallel_RA_inc.h"
#include <time.h>


static void all_to_allv_test(all_to_allv_buffer non_uniform_buffer, int *recv_buffer_offset_size, u64 **recv_buffer, MPI_Comm comm, int *outer_hash_buffer_size, double* at, double* bt, double *avt);

int main(int argc, char **argv)
{
    mpi_comm mcomm;
    mcomm.create(argc, argv);
    srand (time(NULL));

    all_to_allv_buffer non_uniform_buffer;
    non_uniform_buffer.nprocs = mcomm.get_nprocs();
    non_uniform_buffer.ra_count = atoi(argv[1]);
    non_uniform_buffer.local_compute_output_size_total = 0;

    if (mcomm.get_rank() == 0)
        std::cout << "nprocs " << non_uniform_buffer.nprocs << " ra count " << non_uniform_buffer.ra_count <<std::endl;

    non_uniform_buffer.local_compute_output_size_flat = new int[non_uniform_buffer.nprocs * non_uniform_buffer.ra_count];
    non_uniform_buffer.cumulative_tuple_process_map = new int[non_uniform_buffer.nprocs];
    memset(non_uniform_buffer.cumulative_tuple_process_map, 0, non_uniform_buffer.nprocs * sizeof(int));
    non_uniform_buffer.local_compute_output = new vector_buffer*[non_uniform_buffer.ra_count];

    for (int r=0; r < non_uniform_buffer.ra_count; r++)
    {
        non_uniform_buffer.local_compute_output[r] = new vector_buffer[non_uniform_buffer.nprocs];
        for (int i=0; i < non_uniform_buffer.nprocs; i++)
        {
            int random = rand() % 10 + 1;
            non_uniform_buffer.local_compute_output_size_flat[i * non_uniform_buffer.ra_count + r] = atoi(argv[2]) + (random*atoi(argv[2]))/100;
            non_uniform_buffer.cumulative_tuple_process_map[i] = non_uniform_buffer.cumulative_tuple_process_map[i] + non_uniform_buffer.local_compute_output_size_flat[i * non_uniform_buffer.ra_count + r];

            u64 val = i;
            for (int t=0; t < non_uniform_buffer.local_compute_output_size_flat[i * non_uniform_buffer.ra_count + r]; t++)
            {
                non_uniform_buffer.local_compute_output[r][i].vector_buffer_append((const unsigned char*)&val, sizeof(u64));
                non_uniform_buffer.local_compute_output_size_total++;
            }
        }
    }

    int outer_hash_buffer_size = 0;
    u64 *recv_buffer;
    int *recv_buffer_offset_size = new int[non_uniform_buffer.ra_count * non_uniform_buffer.nprocs];
    memset(recv_buffer_offset_size, 0, non_uniform_buffer.ra_count * non_uniform_buffer.nprocs * sizeof(int));

    double at, bt, avt;
    double nu_start = MPI_Wtime();
    all_to_allv_test(non_uniform_buffer, recv_buffer_offset_size, &recv_buffer, MPI_COMM_WORLD, &outer_hash_buffer_size, &at, &bt, &avt);
    double nu_end = MPI_Wtime();

    double max_nu_time = 0;
    double total_nu_time = nu_end - nu_start;
    MPI_Allreduce(&total_nu_time, &max_nu_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (total_nu_time == max_nu_time)
        std::cout << "[NU] Number of processes " << mcomm.get_nprocs() << " Relation count " << non_uniform_buffer.ra_count << " NU time " << (nu_end - nu_start) << " [" << at + bt +avt << "] " << at << " " << bt << " " << avt << std::endl;

    //std::cout << "[Non-uniform] Rank " << mcomm.get_rank() << " buffer size " << outer_hash_buffer_size << std::endl;
    for (int i=0; i < outer_hash_buffer_size; i++)
    {
        if ((int)recv_buffer[i] != mcomm.get_rank())
            std::cout << "U Shout!!!!" << std::endl;
    }

    delete[] recv_buffer;
    delete[] non_uniform_buffer.local_compute_output_size_flat;
    delete[] non_uniform_buffer.cumulative_tuple_process_map;
    for (int i=0; i < non_uniform_buffer.ra_count; i++)
        delete[] non_uniform_buffer.local_compute_output[i];
    delete[] non_uniform_buffer.local_compute_output;
    delete[] recv_buffer_offset_size;


    all_to_all_buffer uniform_buffer;
    uniform_buffer.nprocs = mcomm.get_nprocs();
    uniform_buffer.ra_count = atoi(argv[1]);


    uniform_buffer.local_compute_output = new u64[uniform_buffer.ra_count * uniform_buffer.nprocs * (int)ceil(atoi(argv[2]) * 1.1)];
    for (int i=0; i < uniform_buffer.ra_count * uniform_buffer.nprocs * (int)ceil(atoi(argv[2]) * 1.1); i++)
        uniform_buffer.local_compute_output[i] = i / (uniform_buffer.ra_count * (int)ceil(atoi(argv[2]) * 1.1));


    double u_start = MPI_Wtime();
    u64 *cumulative_all_to_allv_buffer = new u64[uniform_buffer.ra_count * uniform_buffer.nprocs * (int)ceil(atoi(argv[2]) * 1.1)];
    MPI_Alltoall(uniform_buffer.local_compute_output, uniform_buffer.ra_count * (int)ceil(atoi(argv[2]) * 1.1), MPI_UNSIGNED_LONG_LONG, cumulative_all_to_allv_buffer, uniform_buffer.ra_count * (int)ceil(atoi(argv[2]) * 1.1), MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
    double u_end = MPI_Wtime();

    double max_u_time = 0;
    double total_u_time = u_end - u_start;
    MPI_Allreduce(&total_u_time, &max_u_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (total_u_time == max_u_time)
        std::cout << "[U] Number of processes " << mcomm.get_nprocs() << " Relation count " << uniform_buffer.ra_count << " U time " << (u_end - u_start) << std::endl;

    //std::cout << "[Uniform] Rank " << mcomm.get_rank() << " buffer size " << uniform_buffer.ra_count * uniform_buffer.nprocs * (int)ceil(atoi(argv[2]) * 1.1) << std::endl;
    for (int i=0; i < uniform_buffer.ra_count * uniform_buffer.nprocs * (int)ceil(atoi(argv[2]) * 1.1); i++)
    {
        if ((int)cumulative_all_to_allv_buffer[i] != mcomm.get_rank())
            std::cout << "U Shout!!!!" << std::endl;
    }

    delete[] cumulative_all_to_allv_buffer;
    delete[] uniform_buffer.local_compute_output;

    //std::cout << "Number of processes " << mcomm.get_nprocs() << " Relation count " << uniform_buffer.ra_count << " NU time " << (nu_end - nu_start) << " U time " << (u_end - u_start) << std::endl;
    mcomm.destroy();
    return 0;
}


static void all_to_allv_test(all_to_allv_buffer non_uniform_buffer, int *recv_buffer_offset_size, u64 **recv_buffer, MPI_Comm comm, int *outer_hash_buffer_size, double* at, double* bt, double *avt)
{
    u32 RA_count = non_uniform_buffer.ra_count;
    int nprocs = non_uniform_buffer.nprocs;
    int rank;
    MPI_Comm_rank(comm, &rank);


    double at_start = MPI_Wtime();
    MPI_Alltoall(non_uniform_buffer.local_compute_output_size_flat, RA_count, MPI_INT, recv_buffer_offset_size, RA_count, MPI_INT, comm);
    double at_end = MPI_Wtime();
    *at = at_end - at_start;


    double bt_start = MPI_Wtime();
    int *send_disp = new int[nprocs];
    int *recv_counts = new int[nprocs];
    int *recv_displacements = new int[nprocs];

    recv_counts[0] = 0;
    send_disp[0] = 0;
    recv_displacements[0] = 0;

    u64* send_buffer = new u64[non_uniform_buffer.local_compute_output_size_total];

    u32 boffset = 0;
    for(int i = 0; i < nprocs; i++)
    {
        recv_counts[i] = 0;

        if (i >= 1)
            send_disp[i] = send_disp[i - 1] + non_uniform_buffer.cumulative_tuple_process_map[i - 1];

        for (u32 r = 0; r < RA_count; r++)
        {
            memcpy(send_buffer + boffset, non_uniform_buffer.local_compute_output[r][i].buffer, non_uniform_buffer.local_compute_output[r][i].size);
            boffset = boffset + (non_uniform_buffer.local_compute_output[r][i].size)/sizeof(u64);
            non_uniform_buffer.local_compute_output[r][i].vector_buffer_free();

            recv_counts[i] = recv_counts[i] + recv_buffer_offset_size[i*RA_count + r];
        }

        if (i >= 1)
            recv_displacements[i] = recv_displacements[i - 1] + recv_counts[i - 1];
        *outer_hash_buffer_size = *outer_hash_buffer_size + recv_counts[i];
    }

    //std::cout << "X " << non_uniform_buffer.local_compute_output_size_total << " Y " << send_disp[nprocs - 1] << " Z " << non_uniform_buffer.cumulative_tuple_process_map[nprocs - 1] << std::endl;

    //assert(non_uniform_buffer.local_compute_output_size_total == send_disp[nprocs - 1] + non_uniform_buffer.cumulative_tuple_process_map[nprocs - 1]);
#if  1
    *recv_buffer = new u64[*outer_hash_buffer_size];
    double bt_end = MPI_Wtime();
    *bt = bt_end - bt_start;


    double atv_start = MPI_Wtime();
    MPI_Alltoallv(send_buffer, non_uniform_buffer.cumulative_tuple_process_map, send_disp, MPI_UNSIGNED_LONG_LONG, *recv_buffer, recv_counts, recv_displacements, MPI_UNSIGNED_LONG_LONG, comm);
    double atv_end = MPI_Wtime();
    *avt = atv_end - atv_start;

    delete[] send_buffer;
    delete[] send_disp;
    delete[] recv_displacements;
    delete[] recv_counts;
#endif
}

