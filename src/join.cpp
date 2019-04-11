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

double j1x, j2;
double c1 = 0, c2 = 0;
double i1 = 0, i2 = 0;

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


void parallel_read_input_relation_from_file_to_local_buffer(const char *file_name, u64** read_buffer, u32 row_count, u32* local_row_count)
{
    u32 read_offset;
    read_offset = ceil((float)row_count / nprocs) * rank;

    if (read_offset > row_count)
    {
        *local_row_count = 0;
        return;
    }

    if (read_offset + ceil((float)row_count / nprocs) > row_count)
      *local_row_count = row_count - read_offset;
    else
      *local_row_count = (u32) ceil((float)row_count / nprocs);

    if (*local_row_count == 0)
        return;

    char data_filename[1024];
    sprintf(data_filename, "%s/data.raw", file_name);
    //int fp = open(data_filename, O_RDONLY);


    size_t bsize = (u64) *local_row_count *  (u64)COL_COUNT * sizeof(u64);
    *read_buffer = new u64[*local_row_count * COL_COUNT];
    std::ifstream myFile (data_filename, std::ios::in | std::ios::binary);
    myFile.seekg (read_offset * COL_COUNT * sizeof(u64));
    myFile.read ((char*)*read_buffer, bsize);
    myFile.close();

    //u32 rb_size = pread(fp, *read_buffer, *local_row_count * COL_COUNT * sizeof(u64), read_offset * COL_COUNT * sizeof(u64));
    //if (rb_size != *local_row_count * COL_COUNT * sizeof(u64))
    //{
    //    std::cout << "Wrong IO: rank: " << rank << " " << rb_size << " " <<  *local_row_count << " " << COL_COUNT << std::endl;
    //    MPI_Abort(MPI_COMM_WORLD, -1);
    //}
    //close(fp);

    return;
}

void buffer_data_to_hash_buffer(u32 local_number_of_rows, u64* input_data, u64** outer_hash_data, u32* outer_hash_buffer_size, MPI_Comm comm)
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
    for (u32 i = 1; i < nprocs; i++)
    {
        prefix_sum_process_size[i] = prefix_sum_process_size[i - 1] + process_size[i - 1];
        process_size[i] = 2 * process_data_vector[i].size();
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
    for (u32 i = 1; i < nprocs; i++)
    {
        prefix_sum_recv_process_size_buffer[i] = prefix_sum_recv_process_size_buffer[i - 1] + recv_process_size_buffer[i - 1];
        *outer_hash_buffer_size = *outer_hash_buffer_size + recv_process_size_buffer[i];
    }

    *outer_hash_data = new u64[*outer_hash_buffer_size];
    MPI_Alltoallv(process_data, process_size, prefix_sum_process_size, MPI_UNSIGNED_LONG_LONG, *outer_hash_data, recv_process_size_buffer, prefix_sum_recv_process_size_buffer, MPI_UNSIGNED_LONG_LONG, comm);

    *outer_hash_buffer_size = *outer_hash_buffer_size / COL_COUNT;

    delete[] process_data;

    return;
}


u64 parallel_join(relation<2>& O, relation<2>& G, relation<2>& T)
{
    j1x = MPI_Wtime();

    tuple<2> t;
    t[0] = -1;
    t[1] = -1;
    tuple<2> selectall(t);
    relation<2> tempT;

    // Send Join output
    /* process_size[j] stores the number of samples to be sent to process with rank j */
    int process_size[nprocs];
    memset(process_size, 0, nprocs * sizeof(int));

    /* vector[i] contains the data that needs to be sent to process i */
    std::vector<tuple<2>> *process_data_vector;
    process_data_vector = new std::vector<tuple<2>>[nprocs];

    u64 tuple_count = 0;
    u64 non_deduplicate_tuple_count = 0;

    tuple<2> dt;
    tuple<2> s;
    s[1] = -1;
    tuple<2> select;
    uint64_t index;
    for (relation<2>::iter dit(G, selectall); dit.more(); dit.advance())
    {
      s[0] = (*dit)[0];
      select = s;

      for (relation<2>::iter git(T, select); git.more(); git.advance())
      {
        dt[0] = (*dit)[1];
        dt[1] = (*git)[1];
        tuple_count++;

        if (tempT.insert(dt) == true)
        {
            index = outer_hash(dt[1])%nprocs;
            process_data_vector[index].push_back(dt);
        }
      }
    }
    j2 = MPI_Wtime();


    c1 = MPI_Wtime();
    int prefix_sum_process_size[nprocs];
    memset(prefix_sum_process_size, 0, nprocs * sizeof(int));

    process_size[0] = 2 * process_data_vector[0].size();
    for(u32 i = 1; i < nprocs; i++)
    {
        prefix_sum_process_size[i] = prefix_sum_process_size[i - 1] + (2 * process_data_vector[i - 1].size());
        process_size[i] = 2 * process_data_vector[i].size();
        non_deduplicate_tuple_count = non_deduplicate_tuple_count + process_data_vector[i].size();
    }

    int process_data_buffer_size = prefix_sum_process_size[nprocs - 1] + process_size[nprocs - 1];

    u64* process_data = 0;
    process_data = new u64[process_data_buffer_size];
    memset(process_data, 0, process_data_buffer_size * sizeof(u64));

    for(u32 i = 0; i < nprocs; i++)
        memcpy(process_data + prefix_sum_process_size[i], &process_data_vector[i][0], process_data_vector[i].size() * sizeof(tuple<2>));

    delete[] process_data_vector;

    /* This step prepares for actual data transfer */
    /* Every process sends to every other process the amount of data it is going to send */

    int recv_process_size_buffer[nprocs];
    memset(recv_process_size_buffer, 0, nprocs * sizeof(int));
    MPI_Alltoall(process_size, 1, MPI_INT, recv_process_size_buffer, 1, MPI_INT, comm);

    int prefix_sum_recv_process_size_buffer[nprocs];
    memset(prefix_sum_recv_process_size_buffer, 0, nprocs * sizeof(int));

    /* Sending data to all processes: What is the buffer size to allocate */
    u32 outer_hash_buffer_size = recv_process_size_buffer[0];
    for(u32 i = 1; i < nprocs; i++)
    {
        prefix_sum_recv_process_size_buffer[i] = prefix_sum_recv_process_size_buffer[i - 1] + recv_process_size_buffer[i - 1];
        outer_hash_buffer_size = outer_hash_buffer_size + recv_process_size_buffer[i];
    }

    u64 *hash_buffer = 0;
    hash_buffer = new u64[outer_hash_buffer_size];
    memset(hash_buffer, 0, outer_hash_buffer_size * sizeof(u64));

    MPI_Alltoallv(process_data, process_size, prefix_sum_process_size, MPI_UNSIGNED_LONG_LONG, hash_buffer, recv_process_size_buffer, prefix_sum_recv_process_size_buffer, MPI_UNSIGNED_LONG_LONG, comm);


    c2 = MPI_Wtime();

    i1 = MPI_Wtime();
    u64 tcount = 0;

    for (u32 k = 0; k < outer_hash_buffer_size; k = k + 2)
    {
        tuple<2> dt;
        dt[0] = hash_buffer[k];
        dt[1] = hash_buffer[k + 1];

        if (O.insert(dt) == true)
            tcount++;
    }

    delete[] hash_buffer;
    delete[] process_data;
    i2 = MPI_Wtime();


    //if (rank == 0)
    //    std::cout << "Join: " << (j2 - j1x)
    //              << " Comm: " << (c2 - c1)
    //              << " Insert: " << (i2 - i1)
    //              << " Total: " << (i2 - j1x)
    //              << " [" << (j2 - j1x) + (c2 - c1) + (i2 - i1) << "]"
    //              << " Delta: " << tcount
    //              << std::endl;

    return tcount;
}


int main(int argc, char **argv)
{
    // Initializing MPI
    MPI_Init(&argc, &argv);
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    comm = MPI_COMM_WORLD;


    u32 G_global_row_count;
    u32 T_global_row_count;

    nprocs = (u32)size;
    G_global_row_count = atoi(argv[3]);
    T_global_row_count = atoi(argv[4]);


    double G_ior_start = MPI_Wtime();
    u32 G_entry_count;
    u64 *G_input_buffer = NULL;
    parallel_read_input_relation_from_file_to_local_buffer(argv[1], &G_input_buffer, G_global_row_count, &G_entry_count);
    double G_ior_end = MPI_Wtime();

    double G_hash_start = MPI_Wtime();
    u32 G_hash_entry_count;
    u64 *G_hashed_data = NULL;
    buffer_data_to_hash_buffer(G_entry_count, G_input_buffer, &G_hashed_data, &G_hash_entry_count, MPI_COMM_WORLD);
    double G_hash_end = MPI_Wtime();


    double T_ior_start = MPI_Wtime();
    u32 T_entry_count;
    u64 *T_input_buffer = NULL;
    parallel_read_input_relation_from_file_to_local_buffer(argv[2], &T_input_buffer, T_global_row_count, &T_entry_count);
    double T_ior_end = MPI_Wtime();

    double T_hash_start = MPI_Wtime();
    u32 T_hash_entry_count;
    u64 *T_hashed_data = NULL;
    buffer_data_to_hash_buffer(T_entry_count, T_input_buffer, &T_hashed_data, &T_hash_entry_count, MPI_COMM_WORLD);
    double T_hash_end = MPI_Wtime();


    delete[] G_input_buffer;
    delete[] T_input_buffer;


    double T_hash_init_start = MPI_Wtime();
    relation<2> T(COL_COUNT * T_hash_entry_count, T_hashed_data);
    double T_hash_init_end = MPI_Wtime();

    double G_hash_init_start = MPI_Wtime();
    relation<2> G(COL_COUNT * G_hash_entry_count, G_hashed_data);
    double G_hash_init_end = MPI_Wtime();


    delete[] T_hashed_data;
    delete[] G_hashed_data;

    relation<2> O;


    double join_start = MPI_Wtime();
    u64 jcount = parallel_join(O, G, T);

#if 0
    tuple<2> t;
    t[0] = -1;
    t[1] = -1;
    tuple<2> selectall(t);

    for (relation<2>::iter dit(G, selectall); dit.more(); dit.advance())
        std::cout << (*dit)[0] << "\t" << (*dit)[1] << std::endl;

    std::cout << std::endl;

    for (relation<2>::iter dit(T, selectall); dit.more(); dit.advance())
        std::cout << (*dit)[0] << "\t" << (*dit)[1] << std::endl;

    std::cout << std::endl;

    for (relation<2>::iter dit(O, selectall); dit.more(); dit.advance())
        std::cout << (*dit)[0] << "\t" << (*dit)[1] << std::endl;
#endif

    double join_end = MPI_Wtime();

    u64 total_sum = 0;
    MPI_Allreduce(&jcount, &total_sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);

    double elapsed_time = join_end - G_ior_start;
    double max_time = 0;

    MPI_Allreduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (elapsed_time == max_time)
    {
        std::cout << "n: "
                  << nprocs
                  << " rank: "
                  << rank
                  << " G local: "
                  << G_hash_entry_count
                  << " T local: "
                  << T_hash_entry_count
                  << std::endl
                  << "G global: "
                  << G_global_row_count
                  << " T global: "
                  << T_global_row_count
                  << " Join count: " << total_sum
                  << std::endl
                  << "G Read time: " << (G_ior_end - G_ior_start)
                  << " T Read time: " << (T_ior_end - T_ior_start)
                  << " G hash comm: " << (G_hash_end - G_hash_start)
                  << " T hash comm: " << (T_hash_end - T_hash_start)
                  << " G hash init: " << (G_hash_init_end - G_hash_init_start)
                  << " T hash init: " << (T_hash_init_end - T_hash_init_start)
                  << " Join: " << (join_end - join_start)
                  << " Total: " << (join_end - G_ior_start)
                  << " [" << (G_ior_end - G_ior_start) + (T_ior_end - T_ior_start) + (G_hash_end - G_hash_start) + (T_hash_end - T_hash_start) + (G_hash_init_end - G_hash_init_start) + (T_hash_init_end - T_hash_init_start) + (join_end - join_start) << "]" << std::endl;


        std::cout << "Local Join: " << (j2 - j1x)
                  << " Comm: " << (c2 - c1)
                  << " Insert: " << (i2 - i1)
                  << " Total: " << (i2 - j1x)
                  << " [" << (j2 - j1x) + (c2 - c1) + (i2 - i1) << "]"
                  << std::endl;
    }

    // Finalizing MPI
    MPI_Finalize();

    return 0;
}







