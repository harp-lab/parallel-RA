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
#include <unordered_map>
#include <unordered_set>

//#include <omp.h>

#include "btree.h"
#include "btree_relation.h"

static int rank = 0;
static u32 nprocs = 1;
static MPI_Comm comm;
static u32 global_row_count;
static u32 global_col_count;

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

    return h0 ^ h1;// ^ (h1 << 31);
}



u64 outer_hash(const u64 val)
{
    return tunedhash((u8*)(&val),sizeof(u64));
}
#endif


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

    //u64 process_data[process_data_buffer_size];
    u64* process_data = new u64[process_data_buffer_size];

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

    delete[] process_data;
    /*
    if (rank == 0)
        for (u32 i = 0; i < *outer_hash_buffer_size; i=i+2)
            std::cout << "i " << i << " : " << (*outer_hash_data)[i] << ", " << (*outer_hash_data)[i + 1] << std::endl;
            */

    return;
}


std::unordered_map<u64, std::unordered_set<u64>> * parallel_map_join(std::unordered_map<u64, std::unordered_set<u64>>* delT, std::unordered_map<u64, std::unordered_set<u64>>& G, std::unordered_map<u64, std::unordered_set<u64>>& T, int lc, int* lb, int* running_t_count, double* running_time)
{
    std::unordered_map<u64, std::unordered_set<u64>> * delTT = new std::unordered_map<u64, std::unordered_set<u64>>;

    double j1, j2;
    double c1 = 0, c2 = 0;
    double i1 = 0, i2 = 0;
    double v1 = 0;
    double v2 = 0;


    j1 = MPI_Wtime();


    tuple<2> t;
    t[0] = -1;
    t[1] = -1;
    //tuple<2> selectall(t);
    std::unordered_map<u64, std::unordered_set<u64>> tempT;


    // Send Join output
    // process_size[j] stores the number of samples to be sent to process with rank j
    int process_size[nprocs];
    memset(process_size, 0, nprocs * sizeof(int));

    // vector[i] contains the data that needs to be sent to process i
    std::vector<tuple<2>> *process_data_vector;
    process_data_vector = new std::vector<tuple<2>>[nprocs];

    u64 tuple_count = 0;
    u64 non_deduplicate_tuple_count = 0;

    int tcount = 0;

    //tuple<2> dt;
    tuple<2> s;
    s[1] = -1;
    //tuple<2> select;
    uint64_t index;

    tuple<2> dt;
    //u64 k;
#if 0
    std::cout << "----G-------" << std::endl;
    for ( auto local_it = G.begin(); local_it!= G.end(); ++local_it )
    {
        std::unordered_set<u64> k = local_it->second;
        for (auto it2 = k.begin(); it2 != k.end(); it2++)
        {
            std::cout << local_it->first << "\t" << *it2 << std::endl;
        }
    }

    std::cout << "----delT-------" << std::endl;
    for ( auto local_it = delT->begin(); local_it!= delT->end(); ++local_it )
    {
        std::unordered_set<u64> k = local_it->second;
        for (auto it2 = k.begin(); it2 != k.end(); it2++)
        {
            std::cout << local_it->first << "\t" << *it2 << std::endl;
        }
    }
#endif
    for (auto it = delT->begin(); it != delT->end(); it++)
    {
        std::unordered_set<u64> it2 = it->second;
        for (auto dit2 = it2.begin(); dit2 != it2.end(); dit2++)
        {
            std::unordered_set<u64> Git = G[*dit2];
            for (auto it2 = Git.begin(); it2 != Git.end(); it2++)
            {
                tuple_count++;

                auto itx = tempT.find(it->first);
                if( itx != tempT.end() ) {
                    auto it2x = (itx->second).find(*it2);
                    if( it2x != (itx->second).end() ) {
                        ;
                    }
                    else{
                        (itx->second).insert(*it2);
                        tempT[it->first] = itx->second;
                        index = outer_hash(*it2)%nprocs;
                        dt[0] = it->first;
                        dt[1] = *it2;
                        process_data_vector[index].push_back(dt);
                        //row_count++;
                    }
                }
                else {
                    std::unordered_set<u64> k;
                    k.insert(*it2);
                    tempT.insert(std::make_pair(it->first,k));
                    index = outer_hash(*it2)%nprocs;
                    dt[0] = it->first;
                    dt[1] = *it2;
                    process_data_vector[index].push_back(dt);
                    //row_count++;
                }
            }
        }
    }
    j2 = MPI_Wtime();
#if 1


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

#if 1
    MPI_Alltoallv(process_data, process_size, prefix_sum_process_size, MPI_UNSIGNED_LONG_LONG, hash_buffer, recv_process_size_buffer, prefix_sum_recv_process_size_buffer, MPI_UNSIGNED_LONG_LONG, comm);
#endif

#if 0
    MPI_Request* request = new MPI_Request[nprocs * 2];
    MPI_Status* status = new MPI_Status[nprocs * 2];

    for (u32 i = 0; i < nprocs; i++)
    {
        MPI_Irecv(hash_buffer + prefix_sum_recv_process_size_buffer[i], recv_process_size_buffer[i], MPI_UNSIGNED_LONG_LONG, i, 123, MPI_COMM_WORLD, &request[nprocs + i]);
    }

    for (u32 i = 0; i < nprocs; i++)
    {
        MPI_Isend(process_data + prefix_sum_process_size[i], process_size[i], MPI_UNSIGNED_LONG_LONG, i, 123, MPI_COMM_WORLD, &request[i]);
    }

    MPI_Waitall(2*nprocs, request, status);

    delete [] request;
    delete [] status;
#endif


    c2 = MPI_Wtime();
#if 1
    i1 = MPI_Wtime();
    int count = 0;
    u64 tduplicates = 0;
    for(u32 ko = 0; ko < outer_hash_buffer_size; ko = ko + 2)
    {
        auto it = T.find(hash_buffer[ko]);
        if( it != T.end() ) {
            auto it2 = (it->second).find(hash_buffer[ko + 1]);
            if( it2 != (it->second).end() ) {
                tduplicates++;
            }
            else{
                (it->second).insert(hash_buffer[ko + 1]);
                T[hash_buffer[ko]] = it->second;
                tcount++;

                auto sit = (*delTT).find(hash_buffer[ko]);
                if( sit != (*delTT).end() ) {
                    (sit->second).insert(hash_buffer[ko + 1]);
                    (*delTT)[hash_buffer[ko]] = sit->second;
                    count++;
                }
                else {
                    std::unordered_set<u64> k;
                    k.insert(hash_buffer[ko + 1]);
                    (*delTT).insert(std::make_pair(hash_buffer[ko],k));
                    count++;
                }
            }
        }
        else {
            std::unordered_set<u64> k;
            k.insert(hash_buffer[ko + 1]);
            T.insert(std::make_pair(hash_buffer[ko],k));
            tcount++;

            auto sit = (*delTT).find(hash_buffer[ko]);
            if( sit != (*delTT).end() ) {
                (sit->second).insert(hash_buffer[ko + 1]);
                (*delTT)[hash_buffer[ko]] = sit->second;
                count++;
            }
            else {
                std::unordered_set<u64> k;
                k.insert(hash_buffer[ko + 1]);
                (*delTT).insert(std::make_pair(hash_buffer[ko],k));
                count++;
            }
        }
    }

    delete[] hash_buffer;
    delete[] process_data;
    i2 = MPI_Wtime();

    v1 = MPI_Wtime();
    if (lc % 10 == 0)
    {
        int sum = 0;
        MPI_Allreduce(&count, &sum, 1, MPI_INT, MPI_BOR, comm);
        if(sum == 0)
            *lb = 1;
        else
            *lb = 0;
    }
#endif


#endif
    *running_t_count = *running_t_count + tcount;
    v2 = MPI_Wtime();

    *running_time = *running_time + (v2 - j1);

    if (rank == 0)
        std::cout << lc << " [" << v2-j1 << "] [" << *running_time << "]"
                  << " new tuple " << tuple_count
                  << " non duplicate tuples " << non_deduplicate_tuple_count
                  << " Join: " << (j2 - j1)
                  << " Comm: " << (c2 - c1)
                  << " Insert: " << (i2 - i1)
                  << " Verify: " << (v2 - v1)
                  << " Delta: " << tcount
                  << " Duplicte Delta: " << tduplicates
                  << " T : " << *running_t_count
                  << std::endl;

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
    double ior_start = MPI_Wtime();

    nprocs = (u32)size;
    u32 col_count;


    u32 entry_count;
    u64 *input_buffer = NULL;
    parallel_read_input_relation_from_file_to_local_buffer(argv[1], &input_buffer, &entry_count, &col_count);

    double ior_end = MPI_Wtime();

    double hash_start = MPI_Wtime();

    u32 G_hash_entry_count;
    u64 *G_hashed_data = NULL;
    buffer_data_to_hash_buffer(entry_count, col_count, input_buffer,  0 /*hash_column_index*/, &G_hashed_data, &G_hash_entry_count, MPI_COMM_WORLD);

    u32 T_hash_entry_count;
    u64 *T_hashed_data = NULL;
    buffer_data_to_hash_buffer(entry_count, col_count, input_buffer,  1 /*hash_column_index*/, &T_hashed_data, &T_hash_entry_count, MPI_COMM_WORLD);

    delete[] input_buffer;

    double hash_end = MPI_Wtime();


    double relation_start = MPI_Wtime();

    std::unordered_map<u64, std::unordered_set<u64>> T;
    std::unordered_map<u64, std::unordered_set<u64>> G;

    for (u32 i = 0; i < col_count * T_hash_entry_count; i = i + 2)
    {
        auto it = T.find(T_hashed_data[i]);
        if( it != T.end() ) {
            auto it2 = (it->second).find(T_hashed_data[i + 1]);
            if( it2 != (it->second).end() ) {
                ;
            }
            else{
                (it->second).insert(T_hashed_data[i + 1]);
                T[T_hashed_data[i]] = it->second;
            }
        }
        else {
            std::unordered_set<u64> k;
            k.insert(T_hashed_data[i + 1]);
            T.insert(std::make_pair(T_hashed_data[i],k));
        }
    }

    for (u32 i = 0; i < col_count * G_hash_entry_count; i = i + 2)
    {
        auto it = G.find(G_hashed_data[i]);
        if( it != G.end() ) {
            auto it2 = (it->second).find(G_hashed_data[i + 1]);
            if( it2 != (it->second).end() ) {
                ;
            }
            else{
                (it->second).insert(G_hashed_data[i + 1]);
                G[G_hashed_data[i]] = it->second;
            }
        }
        else {
            std::unordered_set<u64> k;
            k.insert(G_hashed_data[i + 1]);
            G.insert(std::make_pair(G_hashed_data[i],k));
        }
    }

    std::unordered_map<u64, std::unordered_set<u64>> * dT = new std::unordered_map<u64, std::unordered_set<u64>>;
    int running_t_count = 0;

    for (u32 i = 0; i < col_count * T_hash_entry_count; i = i + 2)
    {
        auto it = (*dT).find(T_hashed_data[i]);
        if( it != (*dT).end() ) {
            auto it2 = (it->second).find(T_hashed_data[i + 1]);
            if( it2 != (it->second).end() ) {
                ;
            }
            else{
                (it->second).insert(T_hashed_data[i + 1]);
                (*dT)[T_hashed_data[i]] = it->second;
                running_t_count++;
            }
        }
        else {
            std::unordered_set<u64> k;
            k.insert(T_hashed_data[i + 1]);
            (*dT).insert(std::make_pair(T_hashed_data[i],k));
            running_t_count++;
        }
    }

    delete[] T_hashed_data;
    delete[] G_hashed_data;

    std::cout << "[" << rank << "] Initial T count " << running_t_count << std::endl;
    double relation_end = MPI_Wtime();

    double join_start = MPI_Wtime();


    int lb = 0;
    double time = 0;
    dT = parallel_map_join(dT, G, T, 0, &lb, &running_t_count, &time);

    int lc = 1;
    while(true)
    {
      dT = parallel_map_join(dT, G, T, lc, &lb, &running_t_count, &time);

      if (lb == 1)  break;
      lc++;
    }

    delete dT;

    u64 total_sum = 0;
    u64 Tcounter = 0;
    for ( auto local_it = T.begin(); local_it!= T.end(); ++local_it )
    {
        std::unordered_set<u64> k = local_it->second;
        for (auto it2 = k.begin(); it2 != k.end(); it2++)
            Tcounter++;
    }
    MPI_Allreduce(&Tcounter, &total_sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
    double join_end = MPI_Wtime();

    if (rank == 0)
    {
        std::cout << "n: "
                  << nprocs
                  << " G: "
                  << global_row_count
                  << " T: " << total_sum
                  << " Read time: " << (ior_end - ior_start)
                  << " Init hash: " << (hash_end - hash_start)
                  << " Relation init: " << (relation_end - relation_start)
                  << " Fixed point: " << (join_end - join_start)
                  << " Total: " << (join_end - ior_start)
                  << " [" << (ior_end - ior_start) + (hash_end - hash_start) + (relation_end - relation_start) + (join_end - join_start) << "]" << std::endl;
    }

#if 0
    double iow_start = MPI_Wtime();

    char TDname[1024];
    sprintf(TDname, "%s_%d_TC", argv[1], nprocs);
    //std::cout << "Filename " << TCname << std::endl;
    if (rank == 0)
        mkdir(TDname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    MPI_Barrier(MPI_COMM_WORLD);

    char TCname[1024];
    sprintf(TCname, "%s/%s_%d", TDname, argv[1], rank);


    u64 counter = 0;
    u64* buffer = new u64[Tcounter * 2];

    for ( auto local_it = T.begin(); local_it!= T.end(); ++local_it )
    {
        std::unordered_set<u64> k = local_it->second;
        for (auto it2 = k.begin(); it2 != k.end(); it2++)
        {
            memcpy(buffer + counter, &(local_it->first), sizeof(u64));
            counter++;

            memcpy(buffer + counter, &(*it2), sizeof(u64));
            counter++;
        }
    }

    FILE * pFile;
    pFile = fopen (TCname,"wb");
    if (pFile!=NULL)
        fwrite (buffer , sizeof(u64), Tcounter * 2, pFile);
    fclose(pFile);


    //myfile << buffer;
    //myfile.close();
    delete[] buffer;
    double iow_end = MPI_Wtime();


    if (rank == 0)
    {
        std::cout << "n: " << nprocs
                  << " G: " << global_row_count
                  << " T: " << total_sum
                  << " Read time: " << (ior_end - ior_start)
                  << " Init hash: " << (hash_end - hash_start)
                  << " Relation init: " << (relation_end - relation_start)
                  << " Fixed point: " << (join_end - join_start)
                  << " Write time: " << (iow_end - iow_start)
                  << " Total: " << (iow_end - ior_start)
                  << " [" << (ior_end - ior_start) + (hash_end - hash_start) + (relation_end - relation_start) + (join_end - join_start) + (iow_end - iow_start) << "]" << std::endl;
    }
#endif

    // Finalizing MPI
    MPI_Finalize();

    return 0;
}








