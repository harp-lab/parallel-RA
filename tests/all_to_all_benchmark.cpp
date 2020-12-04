#include "../src/parallel_RA_inc.h"
#include <time.h>


static void uniform_benchmark(int ra_count, int nprocs, int epoch_count, int entry_count, int iteration_count);
static void non_uniform_benchmark(int ra_count, int nprocs, int entry_count, int iteration_count);
static void all_to_allv_test(all_to_allv_buffer non_uniform_buffer, MPI_Comm comm, double* at, double* bt, double *avt);

int main(int argc, char **argv)
{
    mpi_comm mcomm;
    mcomm.create(argc, argv);

    srand (time(NULL));


    for (u32 j=4096; j <= 16384; j=j*2)
    {
        for (u32 i= 1; i <= 32; i=i*2)
        {
            if (mcomm.get_rank() == 0)
            {
                std::cout << std::endl;
                std::cout << "[NU] nprocs " << mcomm.get_nprocs() << " ra count " << i << " Entry count " << j <<std::endl;
            }

            non_uniform_benchmark(i, mcomm.get_nprocs(), j, 8);

            for (u32 k= 1; k <= 8; k=k*2)
            {
                if (mcomm.get_rank() == 0)
                {
                    std::cout << std::endl;
                    std::cout << "[U] nprocs " << mcomm.get_nprocs() << " ra count " << i << " Entry count " << j << " Epoch counts " << k << std::endl;
                }
                uniform_benchmark(i, mcomm.get_nprocs(), k, j, 8);
            }
        }
    }


    mcomm.destroy();
    return 0;
}

#if 1
static void uniform_benchmark(int ra_count, int nprocs, int epoch_count, int entry_count, int iteration_count)
{
    all_to_all_buffer uniform_buffer;
    uniform_buffer.nprocs = nprocs;//mcomm.get_nprocs();
    uniform_buffer.ra_count = ra_count;//atoi(argv[1]);

    //std::cout << "BSIZE " <<  (uniform_buffer.ra_count * uniform_buffer.nprocs * (int)ceil(entry_count))/epoch_count << std::endl;
    uniform_buffer.local_compute_output = new u64[(uniform_buffer.ra_count * uniform_buffer.nprocs * (int)ceil(entry_count))/epoch_count];


    u64 *cumulative_all_to_allv_buffer = new u64[(uniform_buffer.ra_count * uniform_buffer.nprocs * (int)ceil(entry_count))/epoch_count];

    for (int it=0; it < iteration_count; it++)
    {
        double u_iter_total=0;
        double u_iter_time[epoch_count];
        double u_start = MPI_Wtime();
        for (int i=0; i<epoch_count; i++)
        {
            double t1 = MPI_Wtime();
            for (int i=0; i < (uniform_buffer.ra_count * uniform_buffer.nprocs * (int)ceil(entry_count)); i=i+epoch_count)
                uniform_buffer.local_compute_output[i/epoch_count] = i / (uniform_buffer.ra_count * (int)ceil(entry_count));
            MPI_Alltoall(uniform_buffer.local_compute_output, (uniform_buffer.ra_count * (int)ceil(entry_count))/epoch_count, MPI_UNSIGNED_LONG_LONG, cumulative_all_to_allv_buffer, (uniform_buffer.ra_count * (int)ceil(entry_count))/epoch_count, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
            double t2 = MPI_Wtime();
            u_iter_time[i] = t2-t1;
            u_iter_total = u_iter_total + (t2-t1);
        }
        double u_end = MPI_Wtime();

        double max_u_time = 0;
        double total_u_time = u_end - u_start;
        MPI_Allreduce(&total_u_time, &max_u_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (total_u_time == max_u_time)
        {
            std::cout << "[U] " << it << " U time " << (u_end - u_start) << " [" << u_iter_total << "] ( ";// << std::endl;
            for (int i=0; i<epoch_count; i++)
                std::cout << u_iter_time[i] << " ";
            std::cout << ")" << std::endl;
        }
    }

    delete[] cumulative_all_to_allv_buffer;
    delete[] uniform_buffer.local_compute_output;
}
#endif


static void non_uniform_benchmark(int ra_count, int nprocs, int entry_count, int iteration_count)
{
    all_to_allv_buffer non_uniform_buffer;
    non_uniform_buffer.nprocs = nprocs;//mcomm.get_nprocs();
    non_uniform_buffer.ra_count = ra_count; //atoi(argv[1]);
    non_uniform_buffer.local_compute_output_size_total = 0;

    non_uniform_buffer.local_compute_output_size_flat = new int[non_uniform_buffer.nprocs * non_uniform_buffer.ra_count];
    non_uniform_buffer.cumulative_tuple_process_map = new int[non_uniform_buffer.nprocs];
    memset(non_uniform_buffer.cumulative_tuple_process_map, 0, non_uniform_buffer.nprocs * sizeof(int));



    for (int it=0; it < iteration_count; it++)
    {
        non_uniform_buffer.local_compute_output = new vector_buffer*[non_uniform_buffer.ra_count];
        for (int r=0; r < non_uniform_buffer.ra_count; r++)
            non_uniform_buffer.local_compute_output[r] = new vector_buffer[non_uniform_buffer.nprocs];

        double at=0, ct=0, bt=0, avt=0, ft=0;

        double nu_start = MPI_Wtime();

        double t1=MPI_Wtime();
        non_uniform_buffer.local_compute_output_size_total=0;
        memset(non_uniform_buffer.cumulative_tuple_process_map, 0, non_uniform_buffer.nprocs * sizeof(int));
        memset(non_uniform_buffer.local_compute_output_size_flat, 0, non_uniform_buffer.nprocs * non_uniform_buffer.ra_count * sizeof(int));
        for (int r=0; r < non_uniform_buffer.ra_count; r++)
        {
            for (int i=0; i < non_uniform_buffer.nprocs; i++)
            {
                int random = rand() % 10 + 1;
                non_uniform_buffer.local_compute_output_size_flat[i * non_uniform_buffer.ra_count + r] = entry_count - (random*entry_count)/100;
                non_uniform_buffer.cumulative_tuple_process_map[i] = non_uniform_buffer.cumulative_tuple_process_map[i] + non_uniform_buffer.local_compute_output_size_flat[i * non_uniform_buffer.ra_count + r];

                u64 val = i;
                for (int t=0; t < non_uniform_buffer.local_compute_output_size_flat[i * non_uniform_buffer.ra_count + r]; t++)
                {
                    non_uniform_buffer.local_compute_output[r][i].vector_buffer_append((const unsigned char*)&val, sizeof(u64));
                    non_uniform_buffer.local_compute_output_size_total++;
                }
            }
        }
        double t2=MPI_Wtime();
        ct = t2-t1;

        all_to_allv_test(non_uniform_buffer, MPI_COMM_WORLD, &at, &bt, &avt);


        t1=MPI_Wtime();
        for (int r=0; r < non_uniform_buffer.ra_count; r++)
            for (int i=0; i < non_uniform_buffer.nprocs; i++)
                non_uniform_buffer.local_compute_output[r][i].vector_buffer_free();
        t2=MPI_Wtime();
        ft = t2-t1;

        double nu_end = MPI_Wtime();

        double max_nu_time = 0;
        double total_nu_time = nu_end - nu_start;
        MPI_Allreduce(&total_nu_time, &max_nu_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (total_nu_time == max_nu_time)
            std::cout << "[NU] " << it << " NU time " << (nu_end - nu_start) << " [" << ct+at + bt +avt +ft << "] " << ct << " " << at << " " << bt << " " << avt << " " << ft << std::endl;

        for (int i=0; i < non_uniform_buffer.ra_count; i++)
            delete[] non_uniform_buffer.local_compute_output[i];
        delete[] non_uniform_buffer.local_compute_output;
    }



    delete[] non_uniform_buffer.local_compute_output_size_flat;
    delete[] non_uniform_buffer.cumulative_tuple_process_map;
}


#if 1
static void all_to_allv_test(all_to_allv_buffer non_uniform_buffer, MPI_Comm comm, double* at, double* bt, double *avt)
{
    double at_start = MPI_Wtime();
    int outer_hash_buffer_size = 0;
    u32 RA_count = non_uniform_buffer.ra_count;
    int nprocs = non_uniform_buffer.nprocs;

    int *recv_buffer_offset_size = new int[RA_count * nprocs];
    memset(recv_buffer_offset_size, 0, RA_count * nprocs * sizeof(int));
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
            recv_counts[i] = recv_counts[i] + recv_buffer_offset_size[i*RA_count + r];
        }

        if (i >= 1)
            recv_displacements[i] = recv_displacements[i - 1] + recv_counts[i - 1];
        outer_hash_buffer_size = outer_hash_buffer_size + recv_counts[i];
    }

    //std::cout << "X " << non_uniform_buffer.local_compute_output_size_total << " Y " << send_disp[nprocs - 1] << " Z " << non_uniform_buffer.cumulative_tuple_process_map[nprocs - 1] << std::endl;

    //assert(non_uniform_buffer.local_compute_output_size_total == send_disp[nprocs - 1] + non_uniform_buffer.cumulative_tuple_process_map[nprocs - 1]);
#if  1
    u64 *recv_buffer = new u64[outer_hash_buffer_size];
    double bt_end = MPI_Wtime();
    *bt = bt_end - bt_start;


    double atv_start = MPI_Wtime();
    MPI_Alltoallv(send_buffer, non_uniform_buffer.cumulative_tuple_process_map, send_disp, MPI_UNSIGNED_LONG_LONG, recv_buffer, recv_counts, recv_displacements, MPI_UNSIGNED_LONG_LONG, comm);


    delete[] recv_buffer;
    delete[] send_buffer;
    delete[] send_disp;
    delete[] recv_displacements;
    delete[] recv_counts;
    delete[] recv_buffer_offset_size;
    double atv_end = MPI_Wtime();
    *avt = atv_end - atv_start;
#endif
}
#endif
