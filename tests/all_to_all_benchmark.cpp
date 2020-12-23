#include "../src/parallel_RA_inc.h"
#include <time.h>
#define ITERATION_COUNT 6

static void uniform_benchmark(int ra_count, int nprocs, int epoch_count, u64 entry_count);
static void non_uniform_benchmark(int ra_count, int nprocs, u64 entry_count, int random_offset, int range);
static void all_to_allv_test(all_to_allv_buffer non_uniform_buffer, MPI_Comm comm, double* all_to_all_time, double* buffer_create_time, double* buffer_delete_time, double *all_to_allv_time, int it, int nprocs, u64 entry_count, int random_offset, int range);

int main(int argc, char **argv)
{
    mpi_comm mcomm;
    mcomm.create(argc, argv);
    srand (time(NULL));
    u32 ra_count = 1;

    for (u64 entry_count=2048; entry_count <= 16384; entry_count=entry_count*2)
        non_uniform_benchmark(ra_count, mcomm.get_nprocs(), entry_count, 90, 10);


    MPI_Barrier(MPI_COMM_WORLD);
    if (mcomm.get_rank() == 0)
        std::cout << "----------------------------------------------------------------" << std::endl << std::endl;

    for (u64 entry_count=2048; entry_count <= 16384; entry_count=entry_count*2)
        non_uniform_benchmark(ra_count, mcomm.get_nprocs(), entry_count, 80, 20);


    MPI_Barrier(MPI_COMM_WORLD);
    if (mcomm.get_rank() == 0)
        std::cout << "----------------------------------------------------------------" << std::endl<< std::endl;

    for (u64 entry_count=2048; entry_count <= 16384; entry_count=entry_count*2)
        non_uniform_benchmark(ra_count, mcomm.get_nprocs(), entry_count, 70, 30);


    MPI_Barrier(MPI_COMM_WORLD);
    if (mcomm.get_rank() == 0)
        std::cout << "----------------------------------------------------------------" << std::endl<< std::endl;

    for (u64 entry_count=2048; entry_count <= 16384; entry_count=entry_count*2)
        non_uniform_benchmark(ra_count, mcomm.get_nprocs(), entry_count, 60, 40);


    MPI_Barrier(MPI_COMM_WORLD);
    if (mcomm.get_rank() == 0)
        std::cout << "----------------------------------------------------------------" << std::endl<< std::endl;

    for (u64 entry_count=2048; entry_count <= 16384; entry_count=entry_count*2)
        non_uniform_benchmark(ra_count, mcomm.get_nprocs(), entry_count, 50, 50);


    MPI_Barrier(MPI_COMM_WORLD);
    if (mcomm.get_rank() == 0)
        std::cout << "----------------------------------------------------------------" << std::endl<< std::endl;

    for (u64 entry_count=2048; entry_count <= 16384; entry_count=entry_count*2)
        non_uniform_benchmark(ra_count, mcomm.get_nprocs(), entry_count, 0, 100);


    MPI_Barrier(MPI_COMM_WORLD);
    if (mcomm.get_rank() == 0)
        std::cout << "----------------------------------------------------------------" << std::endl<< std::endl;

    for (u64 entry_count=2048; entry_count <= 16384; entry_count=entry_count*2)
        for (u32 epoch_count=1; epoch_count<=8; epoch_count=epoch_count*2)
            uniform_benchmark(ra_count, mcomm.get_nprocs(), epoch_count, entry_count);



#if 0
    MPI_Barrier(MPI_COMM_WORLD);

    for (u64 entry_count=4096; entry_count <= 16384; entry_count=entry_count*2)
    {
        for (u32 ra_count= 1; ra_count <= 8; ra_count=ra_count*2)
        {
            if (mcomm.get_rank() == 0)
            {
                std::cout << std::endl;
                std::cout << "---- [NU] nprocs " << mcomm.get_nprocs() << " ra count " << ra_count << " Entry count " << entry_count << " ----" << std::endl;
            }
            non_uniform_benchmark(ra_count, mcomm.get_nprocs(), entry_count);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    for (u64 entry_count=4096; entry_count <= 16384; entry_count=entry_count*2)
    {
        for (u32 rel_count=1; rel_count<=8; rel_count=rel_count*2)
        {
            for (u32 epoch_count=1; epoch_count<=8; epoch_count=epoch_count*2)
            {
                if (mcomm.get_rank() == 0)
                {
                    std::cout << std::endl;
                    std::cout << "[U] nprocs " << mcomm.get_nprocs() << " ra count " << rel_count << " Entry count " << entry_count << " Epoch counts " << epoch_count << " ----" << std::endl;
                }
                uniform_benchmark(rel_count, mcomm.get_nprocs(), epoch_count, entry_count);
            }
        }
    }
#endif
    mcomm.destroy();
    return 0;
}

#if 1
static void uniform_benchmark(int ra_count, int nprocs, int epoch_count, u64 entry_count)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    all_to_all_buffer uniform_buffer;
    uniform_buffer.nprocs = nprocs;//mcomm.get_nprocs();
    uniform_buffer.ra_count = ra_count;//atoi(argv[1]);

    uniform_buffer.local_compute_output = new u64[(uniform_buffer.ra_count * uniform_buffer.nprocs * entry_count)/epoch_count];


    u64 *cumulative_all_to_allv_buffer = new u64[(uniform_buffer.ra_count * uniform_buffer.nprocs * entry_count)/epoch_count];

    for (int it=0; it < ITERATION_COUNT; it++)
    {
        double u_iter_buffer_total=0;
        double u_iter_a2a_total=0;
        double u_iter_buffer_time[epoch_count];
        double u_iter_a2a_time[epoch_count];
        double u_start = MPI_Wtime();

        u64 global_cumulativ_send_count = 0;
        for (int e=0; e<epoch_count; e++)
        {
            double buff_pop_start = MPI_Wtime();
            for (u64 i=0; i < (uniform_buffer.ra_count * uniform_buffer.nprocs * entry_count); i=i+epoch_count)
                uniform_buffer.local_compute_output[i/epoch_count] = i / (uniform_buffer.ra_count * entry_count);
            double buff_pop_end = MPI_Wtime();

            double a2a_start = MPI_Wtime();
            MPI_Alltoall(uniform_buffer.local_compute_output, (uniform_buffer.ra_count * entry_count)/epoch_count, MPI_UNSIGNED_LONG_LONG, cumulative_all_to_allv_buffer, (uniform_buffer.ra_count * entry_count)/epoch_count, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
            double a2a_end = MPI_Wtime();

            if (it == 0)
            {
                u64 global_send_count = 0;
                u64 send_count = ((uniform_buffer.ra_count * entry_count)/epoch_count) * nprocs;
                MPI_Allreduce(&send_count, &global_send_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
                global_cumulativ_send_count = global_cumulativ_send_count + global_send_count;

                if (e == epoch_count - 1)
                {
                    if (rank == 0)
                        std::cout << "[Uniform] Send and Recieve count " << nprocs << " " << epoch_count << " " << entry_count << " " << global_cumulativ_send_count << std::endl;
                    if (global_cumulativ_send_count != entry_count * nprocs * ra_count * nprocs)
                        MPI_Abort(MPI_COMM_WORLD, -1);
                }
            }

            u_iter_buffer_time[e] = buff_pop_end-buff_pop_start;
            u_iter_a2a_time[e] = a2a_end-a2a_start;

            u_iter_buffer_total = u_iter_buffer_total + u_iter_buffer_time[e];
            u_iter_a2a_total = u_iter_a2a_total + u_iter_a2a_time[e];
        }
        double u_end = MPI_Wtime();

        double max_u_time = 0;
        double total_u_time = u_end - u_start;
        MPI_Allreduce(&total_u_time, &max_u_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (total_u_time == max_u_time)
        {
            if (it == 0)
            {
                std::cout << "SKIP [U] " << it << ", " << nprocs << " [ " << entry_count << " " << epoch_count << " ] U time " << (u_end - u_start) << " [" << u_iter_buffer_total + u_iter_a2a_total << "] " << u_iter_buffer_total << " ( ";
            for (int i=0; i<epoch_count; i++)
                std::cout << u_iter_buffer_time[i] << " ";
            std::cout << ") " << u_iter_a2a_total << " ( ";
            for (int i=0; i<epoch_count; i++)
                std::cout << u_iter_a2a_time[i] << " ";
            std::cout << ")" << std::endl;
            }
            else
            {
                std::cout << "[U] " << it << ", " << nprocs << " [ " << entry_count << " " << epoch_count << " ] U time " << (u_end - u_start) << " [" << u_iter_buffer_total + u_iter_a2a_total << "] " << u_iter_buffer_total << " ( ";
            for (int i=0; i<epoch_count; i++)
                std::cout << u_iter_buffer_time[i] << " ";
            std::cout << ") " << u_iter_a2a_total << " ( ";
            for (int i=0; i<epoch_count; i++)
                std::cout << u_iter_a2a_time[i] << " ";
            std::cout << ")" << std::endl;

            }
        }
    }

    delete[] cumulative_all_to_allv_buffer;
    delete[] uniform_buffer.local_compute_output;
}
#endif


static void non_uniform_benchmark(int ra_count, int nprocs, u64 entry_count, int random_offset, int range)
{
    all_to_allv_buffer non_uniform_buffer;
    non_uniform_buffer.nprocs = nprocs;
    non_uniform_buffer.ra_count = ra_count;
    non_uniform_buffer.local_compute_output_size_total = 0;

    non_uniform_buffer.local_compute_output_size_flat = new int[non_uniform_buffer.nprocs * non_uniform_buffer.ra_count];
    non_uniform_buffer.cumulative_tuple_process_map = new int[non_uniform_buffer.nprocs];
    memset(non_uniform_buffer.cumulative_tuple_process_map, 0, non_uniform_buffer.nprocs * sizeof(int));


    for (int it=0; it < ITERATION_COUNT; it++)
    {
        non_uniform_buffer.local_compute_output = new vector_buffer*[non_uniform_buffer.ra_count];
        for (int r=0; r < non_uniform_buffer.ra_count; r++)
            non_uniform_buffer.local_compute_output[r] = new vector_buffer[non_uniform_buffer.nprocs];

        double buffer_creation_time=0, all_to_allv_time=0, buffer_deletion_time=0;

        double nu_start = MPI_Wtime();

        double buffer_creation_start=MPI_Wtime();
        non_uniform_buffer.local_compute_output_size_total=0;
        memset(non_uniform_buffer.cumulative_tuple_process_map, 0, non_uniform_buffer.nprocs * sizeof(int));
        memset(non_uniform_buffer.local_compute_output_size_flat, 0, non_uniform_buffer.nprocs * non_uniform_buffer.ra_count * sizeof(int));

        for (int r=0; r < non_uniform_buffer.ra_count; r++)
        {
            for (int i=0; i < non_uniform_buffer.nprocs; i++)
            {
                int random = random_offset + rand() % range;
                non_uniform_buffer.local_compute_output_size_flat[i * non_uniform_buffer.ra_count + r] = (entry_count * random) / 100;

                non_uniform_buffer.cumulative_tuple_process_map[i] = non_uniform_buffer.cumulative_tuple_process_map[i] + non_uniform_buffer.local_compute_output_size_flat[i * non_uniform_buffer.ra_count + r];
                u64 val = i;
                for (int t=0; t < non_uniform_buffer.local_compute_output_size_flat[i * non_uniform_buffer.ra_count + r]; t++)
                {
                    non_uniform_buffer.local_compute_output[r][i].vector_buffer_append((const unsigned char*)&val, sizeof(u64));
                    non_uniform_buffer.local_compute_output_size_total++;
                }
            }
        }

        double buffer_creation_end=MPI_Wtime();
        buffer_creation_time = buffer_creation_end - buffer_creation_start;

        double all_to_allv_start=MPI_Wtime();
        double a2a_all_to_all_time, a2a_buffer_creation_time, a2a_buffer_deletion_time, a2a_all_to_allv_time;
        all_to_allv_test(non_uniform_buffer, MPI_COMM_WORLD, &a2a_all_to_all_time, &a2a_buffer_creation_time, &a2a_buffer_deletion_time, &a2a_all_to_allv_time, it, nprocs, entry_count, random_offset, range);
        double all_to_allv_end=MPI_Wtime();
        all_to_allv_time = all_to_allv_end - all_to_allv_start;


        double buffer_deletion_start=MPI_Wtime();
        for (int r=0; r < non_uniform_buffer.ra_count; r++)
            for (int i=0; i < non_uniform_buffer.nprocs; i++)
                non_uniform_buffer.local_compute_output[r][i].vector_buffer_free();
        double buffer_deletion_end=MPI_Wtime();
        buffer_deletion_time = buffer_deletion_end-buffer_deletion_start;

        double nu_end = MPI_Wtime();

        double max_nu_time = 0;
        double total_nu_time = nu_end - nu_start;
        MPI_Allreduce(&total_nu_time, &max_nu_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);


        if (max_nu_time == total_nu_time){
            if (it == 0)
                std::cout << "SKIP [NU] " << it << ", " << nprocs << " [ " << random_offset << " " << range << " ] " << entry_count << " NU time " << (nu_end - nu_start) << " [" << buffer_creation_time + all_to_allv_time + buffer_deletion_time << "] Buff pop " << buffer_creation_time << " NA2A " << all_to_allv_time << " [ a2a " << a2a_all_to_all_time << " bct " << a2a_buffer_creation_time << " a2av " << a2a_all_to_allv_time << " bdt " << a2a_buffer_deletion_time << " ] Buff del " << buffer_deletion_time << " " << std::endl;
            else
                std::cout << "[NU] " << it << ", " << nprocs << " [ " << random_offset << " " << range << " ] " << entry_count << " NU time " << (nu_end - nu_start) << " [" << buffer_creation_time + all_to_allv_time + buffer_deletion_time << "] Buff pop " << buffer_creation_time << " NA2A " << all_to_allv_time << " [ a2a " << a2a_all_to_all_time << " bct " << a2a_buffer_creation_time << " a2av " << a2a_all_to_allv_time << " bdt " << a2a_buffer_deletion_time << " ] Buff del " << buffer_deletion_time << " " << std::endl;
        }

        for (int i=0; i < non_uniform_buffer.ra_count; i++)
            delete[] non_uniform_buffer.local_compute_output[i];
        delete[] non_uniform_buffer.local_compute_output;
    }


    delete[] non_uniform_buffer.local_compute_output_size_flat;
    delete[] non_uniform_buffer.cumulative_tuple_process_map;
}


#if 1
static void all_to_allv_test(all_to_allv_buffer non_uniform_buffer, MPI_Comm comm, double* all_to_all_time, double* buffer_create_time, double* buffer_delete_time, double *all_to_allv_time, int it, int nprocs, u64 entry_count, int random_offset, int range)
{
    double all_to_all_start = MPI_Wtime();
    int outer_hash_buffer_size = 0;
    u32 RA_count = non_uniform_buffer.ra_count;
    int *recv_buffer_offset_size = new int[RA_count * nprocs];
    MPI_Alltoall(non_uniform_buffer.local_compute_output_size_flat, RA_count, MPI_INT, recv_buffer_offset_size, RA_count, MPI_INT, comm);
    double all_to_all_end = MPI_Wtime();
    *all_to_all_time = all_to_all_end - all_to_all_start;


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

    u64 *recv_buffer = new u64[outer_hash_buffer_size];
    double bt_end = MPI_Wtime();
    *buffer_create_time = bt_end - bt_start;

    double atv_start = MPI_Wtime();
    MPI_Alltoallv(send_buffer, non_uniform_buffer.cumulative_tuple_process_map, send_disp, MPI_UNSIGNED_LONG_LONG, recv_buffer, recv_counts, recv_displacements, MPI_UNSIGNED_LONG_LONG, comm);
    double atv_end = MPI_Wtime();
    *all_to_allv_time = atv_end - atv_start;

    if (it == 0)
    {
        int rank;
        MPI_Comm_rank(comm, &rank);
        u64 total_send_count=0;
        u64 total_recv_count=0;
        u64 global_total_send_count=0;
        u64 global_total_recv_count=0;
        for (int i=0; i < nprocs; i++)
        {
            total_send_count = total_send_count + non_uniform_buffer.cumulative_tuple_process_map[i];
            total_recv_count = total_recv_count + recv_counts[i];
        }
        MPI_Allreduce(&total_send_count, &global_total_send_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&total_recv_count, &global_total_recv_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        if (rank == 0)
            std::cerr << "[Non Uniform] " << nprocs << " , " << entry_count << " , " << random_offset << " , "  << range << " ,  Total recv count " << global_total_recv_count  << " Total send count " << global_total_send_count << std::endl;

        if (global_total_recv_count != global_total_send_count)
            MPI_Abort(MPI_COMM_WORLD, -1);
        if (global_total_recv_count != global_total_send_count)
            MPI_Abort(MPI_COMM_WORLD, -1);
    }

    double bt_d_start = MPI_Wtime();
    delete[] recv_buffer;
    delete[] send_buffer;
    delete[] send_disp;
    delete[] recv_displacements;
    delete[] recv_counts;
    delete[] recv_buffer_offset_size;
    double bt_d_end = MPI_Wtime();
    *buffer_delete_time = bt_d_end - bt_d_start;
}
#endif
