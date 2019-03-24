#include <chrono>
#include <sys/stat.h>
#include <errno.h>
#include <limits.h>
#include <arpa/inet.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <mpi.h>

#include "btree.h"
#include "btree_relation.h"
//#include "RA.h"

static int rank = 0;
static u32 nprocs = 1;
static MPI_Comm comm;
static u32 global_row_count;
static u32 global_col_count;

uint64_t utime()
{
    return std::chrono::high_resolution_clock::now().time_since_epoch().count();
}


inline static u64 tunedhash(const u8* bp, const u32 len)
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



static u64 outer_hash(const u64 val)
{
    return tunedhash((u8*)(&val),sizeof(u64));
}



void parallel_read_input_relation_from_file_to_local_buffer(const char *file_name, u64** read_buffer, u32* local_row_count, u32* col_count)
{
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

    *col_count = global_col_count;

    int read_offset;
    read_offset = ceil((float)global_row_count / nprocs) * rank;
    if (read_offset + ceil((float)global_row_count / nprocs) > global_row_count)
      *local_row_count = global_row_count - read_offset;
    else
      *local_row_count = (int) ceil((float)global_row_count / nprocs);

    if (*local_row_count < 0)
        *local_row_count = 0;

    char data_filename[1024];
    sprintf(data_filename, "%s/data.raw", file_name);
    int fp = open(data_filename, O_RDONLY);

    //std::cout << "Global row count " << global_row_count << " Col count " << global_col_count << std::endl;
    //std::cout << "Local row count " << *local_row_count << std::endl;

    *read_buffer = new u64[*local_row_count * global_col_count];
    u32 rb_size = pread(fp, *read_buffer, *local_row_count * global_col_count * sizeof(u64), read_offset * global_col_count * sizeof(u64));
    if (rb_size != *local_row_count * global_col_count * sizeof(u64))
    {
        printf("Wrong input format (Meta Data)\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    close(fp);

    //if (rank == 0)
    //printf("Rank %d reads %d elements from %d offset from %s\n", rank, *local_row_count, read_offset, data_filename);

    return;
}

void buffer_data_to_hash_buffer(u32 local_number_of_rows, int col_count, u64* input_data,  int hash_column_index, u64** outer_hash_data, u32* outer_hash_buffer_size, MPI_Comm comm)
{
    /* process_size[j] stores the number of samples to be sent to process with rank j */
    int process_size[nprocs];
    memset(process_size, 0, nprocs * sizeof(int));

    /* vector[i] contains the data that needs to be sent to process i */
    std::vector<u64> process_data_vector[nprocs];

    if (hash_column_index == 0)
    {
        for (u32 i = 0; i < local_number_of_rows * col_count; i=i+2)
        {
            uint64_t index = outer_hash(input_data[i])%nprocs;
            process_size[index] = process_size[index] + col_count;

            process_data_vector[index].push_back(input_data[i]);
            process_data_vector[index].push_back(input_data[i + 1]);

            //if (rank == 0)
            //    std::cout << input_data[i] << " : " << input_data[i + 1] << std::endl;
        }
    }
    else
    {
        for (u32 i = 0; i < local_number_of_rows * col_count; i=i+2)
        {
            uint64_t index = outer_hash(input_data[i+1])%nprocs;
            process_size[index] = process_size[index] + col_count;

            process_data_vector[index].push_back(input_data[i]);
            process_data_vector[index].push_back(input_data[i + 1]);
        }
    }


    int prefix_sum_process_size[nprocs];
    memset(prefix_sum_process_size, 0, nprocs * sizeof(int));
    for (u32 i = 1; i < nprocs; i++)
        prefix_sum_process_size[i] = prefix_sum_process_size[i - 1] + process_size[i - 1];

    /*
    if (rank == 0)
        for (int i = 0; i < nprocs; i++)
            std::cout << "i " << i << " : " << process_size[i]  << std::endl;
    */

    int process_data_buffer_size = prefix_sum_process_size[nprocs - 1] + process_size[nprocs - 1];
    u64 process_data[process_data_buffer_size];

    for (u32 i = 0; i < nprocs; i++)
        memcpy(process_data + prefix_sum_process_size[i], &process_data_vector[i][0], process_data_vector[i].size() * sizeof(u64));

    /*
    if (rank == 0)
        for (int i = 0; i < process_data_buffer_size; i=i+2)
            std::cout << "i " << i << " : " << process_data[i] << ", " << process_data[i + 1] << std::endl;
            */

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
    *outer_hash_buffer_size = 0;
    for(u32 i = 0; i < nprocs; i++)
        *outer_hash_buffer_size = *outer_hash_buffer_size + recv_process_size_buffer[i];

#if 1
    uint total_row_size = 0;
    MPI_Allreduce(outer_hash_buffer_size, &total_row_size, 1, MPI_INT, MPI_SUM, comm);
    if(total_row_size != global_row_count * global_col_count)
    {
        printf("Incorrect distribution %d != %d %d\n", total_row_size, global_row_count, global_col_count);
        MPI_Abort(comm, -1);
    }
#endif

    *outer_hash_data = new u64[*outer_hash_buffer_size];
    MPI_Alltoallv(process_data, process_size, prefix_sum_process_size, MPI_UNSIGNED_LONG_LONG, *outer_hash_data, recv_process_size_buffer, prefix_sum_recv_process_size_buffer, MPI_UNSIGNED_LONG_LONG, comm);

    *outer_hash_buffer_size = *outer_hash_buffer_size / col_count;

    /*
    if (rank == 0)
        for (u32 i = 0; i < *outer_hash_buffer_size; i=i+2)
            std::cout << "i " << i << " : " << (*outer_hash_data)[i] << ", " << (*outer_hash_data)[i + 1] << std::endl;
            */

    return;
}


relation<2> * parallel_join(relation<2>* delT, relation<2>& G, relation<2>& T, int lc, int* lb, int* running_t_count, u64* running_time)
{
    u64 start = utime();
    tuple<2> t;
    t[0] = -1;
    t[1] = -1;
    tuple<2> selectall(t);
    relation<2> tempT;
    relation<2> * delTT = new relation<2>;

    // Send Join output
    /* process_size[j] stores the number of samples to be sent to process with rank j */
    int process_size[nprocs];
    memset(process_size, 0, nprocs * sizeof(int));

    /* vector[i] contains the data that needs to be sent to process i */
    std::vector<u64> *process_data_vector;
    process_data_vector = new std::vector<u64>[nprocs];

    int count = 0;
    int tcount = 0;
    int tuple_count = 0;

    for (relation<2>::iter dit(*delT, selectall); dit.more(); dit.advance())
    {
      tuple<2> s;
      s[0] = (*dit)[1];
      s[1] = -1;
      tuple<2> select(s);

      for (relation<2>::iter git(G, select); git.more(); git.advance())
      {
        tuple<2> dt;
        dt[0] = (*dit)[0];
        dt[1] = (*git)[1];
        tuple_count++;

        if (tempT.insert(dt) == true)
        {
            uint64_t index = outer_hash(dt[1])%nprocs;
            process_size[index] = process_size[index] + 2/*COL_COUNT*/;
            process_data_vector[index].push_back(dt[0]);
            process_data_vector[index].push_back(dt[1]);
        }
      }
    }

    int prefix_sum_process_size[nprocs];
    memset(prefix_sum_process_size, 0, nprocs * sizeof(int));
    for(u32 i = 1; i < nprocs; i++)
        prefix_sum_process_size[i] = prefix_sum_process_size[i - 1] + process_size[i - 1];

    int process_data_buffer_size = prefix_sum_process_size[nprocs - 1] + process_size[nprocs - 1];

    u64* process_data = 0;
    process_data = new u64[process_data_buffer_size];
    memset(process_data, 0, process_data_buffer_size * sizeof(u64));

    for(u32 i = 0; i < nprocs; i++)
        memcpy(process_data + prefix_sum_process_size[i], &process_data_vector[i][0], process_data_vector[i].size() * sizeof(u64));

    delete[] process_data_vector;


    /* This step prepares for actual data transfer */
    /* Every process sends to every other process the amount of data it is going to send */

    int recv_process_size_buffer[nprocs];
    memset(recv_process_size_buffer, 0, nprocs * sizeof(int));
    MPI_Alltoall(process_size, 1, MPI_INT, recv_process_size_buffer, 1, MPI_INT, comm);

    int prefix_sum_recv_process_size_buffer[nprocs];
    memset(prefix_sum_recv_process_size_buffer, 0, nprocs * sizeof(int));
    for(u32 i = 1; i < nprocs; i++)
        prefix_sum_recv_process_size_buffer[i] = prefix_sum_recv_process_size_buffer[i - 1] + recv_process_size_buffer[i - 1];

    /* Sending data to all processes */
    /* What is the buffer size to allocate */
    u32 outer_hash_buffer_size = 0;
    for(u32 i = 0; i < nprocs; i++)
        outer_hash_buffer_size = outer_hash_buffer_size + recv_process_size_buffer[i];

    u64 *hash_buffer = 0;
    hash_buffer = new u64[outer_hash_buffer_size];
    memset(hash_buffer, 0, outer_hash_buffer_size * sizeof(u64));

    MPI_Alltoallv(process_data, process_size, prefix_sum_process_size, MPI_UNSIGNED_LONG_LONG, hash_buffer, recv_process_size_buffer, prefix_sum_recv_process_size_buffer, MPI_UNSIGNED_LONG_LONG, comm);

    count = 0;
    for(u32 k = 0; k < outer_hash_buffer_size; k = k + 2)
    {
        tuple<2> dt;
        dt[0] = hash_buffer[k];
        dt[1] = hash_buffer[k + 1];

        if (T.insert(dt) == true)
        {
            tcount++;
            if (delTT->insert(dt) == true)
              count++;
        }
    }

    int sum = 0;
    MPI_Allreduce(&count, &sum, 1, MPI_INT, MPI_BOR, comm);
    if(sum == 0)
      *lb = 1;
    else
      *lb = 0;

    *running_t_count = *running_t_count + tcount;

    u64 end = utime();
    u64 dTime = (end - start) / 1000000;
    *running_time = *running_time + dTime;

    if (rank == 0)
        std::cout << lc << " [" << dTime << "]  [" << *running_time << "] Tuple count " << tuple_count << " Delta count " << count << " T count: " << *running_t_count << " : " << std::endl;

    delete delT;

    return delTT;
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
    u32 col_count;


    u32 entry_count;
    u64 *input_buffer = NULL;
    parallel_read_input_relation_from_file_to_local_buffer(argv[1], &input_buffer, &entry_count, &col_count);

    u32 G_hash_entry_count;
    u64 *G_hashed_data = NULL;
    buffer_data_to_hash_buffer(entry_count, col_count, input_buffer,  0 /*hash_column_index*/, &G_hashed_data, &G_hash_entry_count, MPI_COMM_WORLD);


    u32 T_hash_entry_count;
    u64 *T_hashed_data = NULL;
    buffer_data_to_hash_buffer(entry_count, col_count, input_buffer,  1 /*hash_column_index*/, &T_hashed_data, &T_hash_entry_count, MPI_COMM_WORLD);

    delete[] input_buffer;

    //std::cout << "[" << rank << "] G Count : " << G_hash_entry_count << " T Count : " << T_hash_entry_count << std::endl;

    /*
    if (rank == 0)
    for (u32 i = 0; i < hash_entry_count * col_count; i=i+2)
        std::cout << hashed_data[i] << " " << hashed_data[i+1] << std::endl;
    */

    relation<2> G(col_count * G_hash_entry_count, G_hashed_data);
    relation<2> T(col_count * T_hash_entry_count, T_hashed_data);

    delete[] T_hashed_data;
    delete[] G_hashed_data;

    relation<2> * dT = new relation<2>;


    tuple<2> t;
    t[0] = -1; t[1] = -1;
    tuple<2> selectall(t);

    int running_t_count = 0;
    for (relation<2>::iter it(T, selectall); it.more(); it.advance())
    {
      tuple<2> t1;
      t1[0] = (*it)[0];
      t1[1] = (*it)[1];
      if (dT->insert(t1) == true)
          running_t_count++;
    }
    //std::cout << "[" << rank << "] Initial T count " << running_t_count << std::endl;

#if 1
    int lc = 0;
    int lb = 0;
    u64 time = 0;

    //std::cout << "Loop count ";
    dT = parallel_join(dT, G, T, 0, &lb, &running_t_count, &time);

    lc++;
    while(true)
    {
      //std::cout << "Loop count ";
      dT = parallel_join(dT, G, T, lc, &lb, &running_t_count, &time);

      if (lb == 1)  break;
      lc++;
    }

    delete dT;

    //char TCname[1024];
    //sprintf(TCname, "%s_TC4", argv[1]);
    //std::cout << "Filename " << TCname << std::endl;

    //std::ofstream myfile;
    //myfile.open (TCname);
    u32 Tcounter = 0;
    for (relation<2>::iter Tit(T, selectall); Tit.more(); Tit.advance())
    {
        Tcounter++;
        //myfile << (*Tit)[0] << "\t" << (*Tit)[1] << "\n";
    }
    //myfile.close();

    u32 total_sum = 0;
    MPI_Allreduce(&Tcounter, &total_sum, 1, MPI_INT, MPI_SUM, comm);

    if (rank == 0)
        std::cout << "Total tuple count: " << total_sum << std::endl;

#endif
    // Finalizing MPI
    MPI_Finalize();

    return 0;
}







