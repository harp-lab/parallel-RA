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
#include "compat.h"
#include "tuple.h"
#include "vector.h"
#include <unordered_set>

//#define DEBUG 1

#define COL_COUNT 2

#include "btree/btree_map.h"

typedef btree::btree_map<u64, u64> Relation0Map;
typedef btree::btree_map<u64, btree::btree_map<u64, u64>* > Relation1Map;

static int rank = 0;
static u32 nprocs = 1;
static MPI_Comm comm;
static u32 global_row_count;



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



u64 hash_function(const u64 val)
{
    return tunedhash((u8*)(&val),sizeof(u64));
}


void parallel_read_input_relation_from_file_to_local_buffer(const char *file_name, u64** read_buffer, u32* local_row_count)
{
    if (rank == 0)
    {
        int global_col_count = 0;
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
        assert(global_col_count == 2);
        fclose(fp_in);
    }

    MPI_Bcast(&global_row_count, 1, MPI_INT, 0, MPI_COMM_WORLD);


    /*
    int read_offset;
    read_offset = ceil((float)global_row_count / nprocs) * rank;
    if (read_offset + ceil((float)global_row_count / nprocs) > global_row_count)
        *local_row_count = global_row_count - read_offset;
    else
        *local_row_count = (int) ceil((float)global_row_count / nprocs);

    if (*local_row_count < 0)
        *local_row_count = 0;
    */

    u32 read_offset;
    read_offset = ceil((float)global_row_count / nprocs) * rank;
    if (read_offset > global_row_count)
    {
        *local_row_count = 0;
        return;
    }
    if (read_offset + ceil((float)global_row_count / nprocs) > global_row_count)
        *local_row_count = global_row_count - read_offset;
    else
        *local_row_count = (u32) ceil((float)global_row_count / nprocs);
    if (*local_row_count == 0)
        return;

    char data_filename[1024];
    sprintf(data_filename, "%s/data.raw", file_name);
    int fp = open(data_filename, O_RDONLY);

    //std::cout << "Global row count " << global_row_count << " Col count " << COL_COUNT << std::endl;
    //std::cout << "Local row count " << *local_row_count << std::endl;

    *read_buffer = new u64[*local_row_count * COL_COUNT];
    u64 rb_size = pread(fp, *read_buffer, *local_row_count * COL_COUNT * sizeof(u64), read_offset * COL_COUNT * sizeof(u64));
    if (rb_size != *local_row_count * COL_COUNT * sizeof(u64))
    {
        printf("Wrong input format (Meta Data)\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    close(fp);

    //if (rank == 0)
    //printf("Rank %d reads %d elements from %d offset from %s\n", rank, *local_row_count, read_offset, data_filename);

    return;
}

void buffer_data_to_hash_buffer(u32 local_number_of_rows, u64* input_data,  int hash_column_index, u64** outer_hash_data, u32* outer_hash_buffer_size, MPI_Comm comm, int buckets, u32* subbuckets_G, u32* subbuckets_T, int** gmap_sub_bucket_rank)
{
    /* process_size[j] stores the number of samples to be sent to process with rank j */
    int process_size[nprocs];
    memset(process_size, 0, nprocs * sizeof(int));

    /* vector[i] contains the data that needs to be sent to process i */
    //std::vector<u64> process_data_vector[nprocs];
    vector_buffer* process_data_vector = (vector_buffer*)malloc(sizeof(vector_buffer) * nprocs);
    for (u32 i = 0; i < nprocs; ++i) {
        process_data_vector[i] = vector_buffer_create_empty();
    }

    u64 val[2] = {0, 0};
    if (hash_column_index == 0)
    {
        for (u32 i = 0; i < local_number_of_rows * COL_COUNT; i=i+2)
        {
            uint64_t bucket_id = hash_function(input_data[i]) % buckets;
            uint64_t sub_bucket_id = hash_function(input_data[i+1]) % subbuckets_G[bucket_id];

            //uint64_t index = (hash_function((bucket_id << 32)^sub_bucket_id))%nprocs;
            int index = gmap_sub_bucket_rank[bucket_id][sub_bucket_id];
            process_size[index] = process_size[index] + COL_COUNT;

            //process_data_vector[index].push_back(input_data[i]);
            //process_data_vector[index].push_back(input_data[i + 1]);

            val[0] = input_data[i];
            val[1] = input_data[i + 1];
            vector_buffer_append(&process_data_vector[index],(unsigned char *) val, sizeof(u64)*2);
        }
    }
    else
    {
        for (u32 i = 0; i < local_number_of_rows * COL_COUNT; i=i+2)
        {
            uint64_t bucket_id = hash_function(input_data[i+1]) % buckets;
            uint64_t sub_bucket_id = hash_function(input_data[i]) % subbuckets_T[bucket_id];

            //uint64_t index = (hash_function((bucket_id << 32)^sub_bucket_id))%nprocs;
            int index = gmap_sub_bucket_rank[bucket_id][sub_bucket_id];
            process_size[index] = process_size[index] + COL_COUNT;

            //process_data_vector[index].push_back(input_data[i]);
            //process_data_vector[index].push_back(input_data[i + 1]);

            val[0] = input_data[i];
            val[1] = input_data[i + 1];
            vector_buffer_append(&process_data_vector[index],(unsigned char *) val, sizeof(u64)*2);
        }
    }

    int prefix_sum_process_size[nprocs];
    memset(prefix_sum_process_size, 0, nprocs * sizeof(int));

    for (u32 i = 1; i < nprocs; i++)
        prefix_sum_process_size[i] = prefix_sum_process_size[i - 1] + process_size[i - 1];

    int process_data_buffer_size = prefix_sum_process_size[nprocs - 1] + process_size[nprocs - 1];
    u64* process_data = new u64[process_data_buffer_size];

    for (u32 i = 0; i < nprocs; i++)
    {
        memcpy(process_data + prefix_sum_process_size[i], (&process_data_vector[i])->buffer, (&process_data_vector[i])->size);
        vector_buffer_free(&process_data_vector[i]);
    }
    free (process_data_vector);

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
    if(total_row_size != global_row_count * COL_COUNT)
    {
        printf("Incorrect distribution %d != %d %d\n", total_row_size, global_row_count, COL_COUNT);
        MPI_Abort(comm, -1);
    }
#endif

    *outer_hash_data = new u64[*outer_hash_buffer_size];
    MPI_Alltoallv(process_data, process_size, prefix_sum_process_size, MPI_UNSIGNED_LONG_LONG, *outer_hash_data, recv_process_size_buffer, prefix_sum_recv_process_size_buffer, MPI_UNSIGNED_LONG_LONG, comm);

    *outer_hash_buffer_size = *outer_hash_buffer_size / COL_COUNT;

    delete[] process_data;

    return;
}


void parallel_map_join(Relation1Map*& delT, Relation1Map*& G, u32* gmap_bucket, Relation1Map*& T, u32* tmap, int lc, int* lb, int* running_t_count, double* running_time, u32 buckets, u32* subbuckets_T, u32* g_t_bucket_indices, u32 g_t_actual_bucket_count, int* gmap_distinct_sub_bucket_rank_count,  int** gmap_distinct_sub_bucket_rank, int* tmap_distinct_sub_bucket_rank_count, int** tmap_distinct_sub_bucket_rank, int** tmap_sub_bucket_rank, u32**& tmap_sub_bucket)
{

    double c = 0;
    double jo = 0;

    double t1_s = MPI_Wtime();
    double c_j = 0;

    u64 val[2] = {0, 0};

    //std::vector<tuple<2>> *process_data_vector;
    //process_data_vector = new std::vector<tuple<2>>[nprocs];
    vector_buffer* process_data_vector = (vector_buffer*)malloc(sizeof(vector_buffer) * nprocs);
    for (u32 i = 0; i < nprocs; ++i) {
        process_data_vector[i] = vector_buffer_create_empty();
    }

    u64 tuple_count = 0;


#if 1
    int process_size[nprocs];
    memset(process_size, 0, nprocs * sizeof(int));
    for (u32 i3 = 0; i3 < g_t_actual_bucket_count; i3++)
    {
        double c_s = MPI_Wtime();
        u32 i = g_t_bucket_indices[i3];

#if DEBUG
        if (rank == 0)
            std::cout << "[" << i << "] Meta Comm Start" << std::endl;
#endif

        //std::vector<tuple<2>> delT_temp;
        vector_buffer delT_temp = vector_buffer_create_empty();

        for (auto it = delT[i].begin(); it != delT[i].end(); it++)
        {
            Relation0Map* it2 = it->second;
            for (auto dit2 = it2->begin(); dit2 != it2->end(); dit2++)
            {
                val[0] = it->first;
                val[1] = dit2->first;

                vector_buffer_append(&delT_temp,(unsigned char *) val, sizeof(u64)*2);
            }
        }

        int buffer_size = (&delT_temp)->size / sizeof(u64);
        u32* size_buf = new u32[tmap_distinct_sub_bucket_rank_count[i]];
        memset(size_buf, 0, sizeof(u32) * tmap_distinct_sub_bucket_rank_count[i]);

        u32 req_counter1 = 0;
        MPI_Request *req1 = new MPI_Request[gmap_distinct_sub_bucket_rank_count[i] + tmap_distinct_sub_bucket_rank_count[i]];
        MPI_Status *stat1 = new MPI_Status[gmap_distinct_sub_bucket_rank_count[i] + tmap_distinct_sub_bucket_rank_count[i]];

        // meta-data
        if (tmap[i] == 1)
        {
            for (int r = 0; r < gmap_distinct_sub_bucket_rank_count[i]; r++)
            {
                MPI_Isend(&buffer_size, 1, MPI_INT, gmap_distinct_sub_bucket_rank[i][r], 123, MPI_COMM_WORLD, &req1[req_counter1]);
                req_counter1++;
            }
        }

        if (gmap_bucket[i] == 1)
        {
            for (int r = 0; r < tmap_distinct_sub_bucket_rank_count[i]; r++)
            {
                MPI_Irecv(size_buf + r, 1, MPI_INT, tmap_distinct_sub_bucket_rank[i][r], 123, MPI_COMM_WORLD, &req1[req_counter1]);
                req_counter1++;
            }
        }

        MPI_Waitall(req_counter1, req1, stat1);

#if DEBUG
        if (rank == 0)
            std::cout << "[" << i << "] Meta Comm End" << std::endl;

        if (rank == 0)
            std::cout << "[" << i << "] Comm Start" << std::endl;
#endif

        u64 total_buffer_size = 0;
        for (int r = 0; r < tmap_distinct_sub_bucket_rank_count[i]; r++)
            total_buffer_size = total_buffer_size + size_buf[r];

#if 1
        u32 req_counter2 = 0;
        MPI_Request *req2 = new MPI_Request[gmap_distinct_sub_bucket_rank_count[i] + tmap_distinct_sub_bucket_rank_count[i]];
        MPI_Status *stat2 = new MPI_Status[gmap_distinct_sub_bucket_rank_count[i] + tmap_distinct_sub_bucket_rank_count[i]];


        // data
        if (tmap[i] == 1)
        {
            for (int r = 0; r < gmap_distinct_sub_bucket_rank_count[i]; r++)
            {
                if (buffer_size != 0)
                {
                    MPI_Isend((&delT_temp)->buffer, buffer_size, MPI_UNSIGNED_LONG_LONG, gmap_distinct_sub_bucket_rank[i][r], 123, MPI_COMM_WORLD, &req2[req_counter2]);
                    req_counter2++;
                }
            }
        }

        u32 offset = 0;
        u64 *recvbuf = new u64[total_buffer_size];
        if (gmap_bucket[i] == 1)
        {
            for (int r = 0; r < tmap_distinct_sub_bucket_rank_count[i]; r++)
            {
                if (size_buf[r] != 0)
                {
                    MPI_Irecv(recvbuf + offset, size_buf[r], MPI_UNSIGNED_LONG_LONG, tmap_distinct_sub_bucket_rank[i][r], 123, MPI_COMM_WORLD, &req2[req_counter2]);
                    offset = offset + size_buf[r];
                    req_counter2++;
                }
            }
        }

        MPI_Waitall(req_counter2, req2, stat2);
        vector_buffer_free(&delT_temp);

        delete[] req1;
        delete[] req2;

        delete[] stat1;
        delete[] stat2;
        delete[] size_buf;

#if DEBUG
        if (rank == 0)
            std::cout << "[" << i << "] Comm End" << std::endl;
#endif
        double c_e = MPI_Wtime();
        c = c + (c_e - c_s);



        double jo_s = MPI_Wtime();

#if DEBUG
        if (rank == 0)
            std::cout << "[" << i << "] Join Start" << std::endl;
#endif

        Relation1Map tempT;
        for (u32 k1 = 0; k1 < total_buffer_size; k1=k1+2)
        {
            auto itd = G[i].find(recvbuf[k1+1]);
            if( itd != G[i].end() ) {
                Relation0Map* Git = itd->second;
                for (auto it2 = Git->begin(); it2 != Git->end(); it2++)
                {
                    tuple_count++;

                    auto itx = tempT.find(recvbuf[k1]);
                    if( itx != tempT.end() ) {
                        auto it2x = (itx->second)->find(it2->first);
                        if( it2x != (itx->second)->end() ) {
                            ;
                        }
                        else{
                            (itx->second)->insert(std::make_pair(it2->first, 0));
                            tempT[recvbuf[k1]] = itx->second;

                            // because we are sending to T and dt and they are hashed on the second column
                            uint64_t bucket_id = hash_function(it2->first) % buckets;
                            uint64_t sub_bucket_id = hash_function(recvbuf[k1]) % subbuckets_T[bucket_id];

                            int index = tmap_sub_bucket_rank[bucket_id][sub_bucket_id];
                            //uint64_t index = (hash_function((bucket_id << 32)^sub_bucket_id))%nprocs;

                            val[0] = recvbuf[k1];
                            val[1] = it2->first;
                            vector_buffer_append(&process_data_vector[index], (unsigned char *) val, sizeof(u64)*2);

                            process_size[index] = process_size[index] + COL_COUNT;
                            //process_data_vector[index].push_back(dt);
                        }
                    }
                    else {
                        Relation0Map* k = new Relation0Map();
                        k->insert(std::make_pair(it2->first, 0));
                        //tempT.insert(std::make_pair(it->first,k));
                        tempT[recvbuf[k1]] = k;

                        // because we are sending to T and dt and they are hashed on the second column
                        uint64_t bucket_id = hash_function(it2->first) % buckets;
                        uint64_t sub_bucket_id = hash_function(recvbuf[k1]) % subbuckets_T[bucket_id];

                        //uint64_t index = (hash_function((bucket_id << 32)^sub_bucket_id))%nprocs;
                        int index = tmap_sub_bucket_rank[bucket_id][sub_bucket_id];

                        val[0] = recvbuf[k1];
                        val[1] = it2->first;
                        vector_buffer_append(&process_data_vector[index], (unsigned char *) val, sizeof(u64)*2);

                        process_size[index] = process_size[index] + COL_COUNT;
                        //process_data_vector[index].push_back(dt);
                    }
                }
            }
        }

        Relation1Map::iterator ix = tempT.begin();
        for(; ix != tempT.end(); ix++)
            delete (ix->second);
        delete[]  recvbuf;

#if DEBUG
        if (rank == 0)
            std::cout << "[" << i << "] Join End" << std::endl;
#endif

        double jo_e = MPI_Wtime();
        jo = jo + (jo_e - jo_s);
#endif
    }
#endif
    double t1_e = MPI_Wtime();
    c_j = (t1_e - t1_s);



#if 1
#if 1
    double c1 = MPI_Wtime();

#if DEBUG
    if (rank == 0)
        std::cout << "Global Comm Start" << std::endl;
#endif

    int prefix_sum_process_size[nprocs];
    memset(prefix_sum_process_size, 0, nprocs * sizeof(int));

    //u64 non_deduplicate_tuple_count = 0;
    //non_deduplicate_tuple_count = non_deduplicate_tuple_count + process_data_vector[0].size();
    for (u32 i = 1; i < nprocs; i++)
    {
        prefix_sum_process_size[i] = prefix_sum_process_size[i - 1] + process_size[i - 1];// (2 * process_data_vector[i - 1].size());
        //non_deduplicate_tuple_count = non_deduplicate_tuple_count + process_data_vector[i].size();
    }

    int process_data_buffer_size = prefix_sum_process_size[nprocs - 1] + process_size[nprocs - 1];

    u64* process_data = 0;
    process_data = new u64[process_data_buffer_size];
    memset(process_data, 0, process_data_buffer_size * sizeof(u64));

    for(u32 i = 0; i < nprocs; i++)
    {
        //memcpy(process_data + prefix_sum_process_size[i], &process_data_vector[i][0], process_data_vector[i].size() * sizeof(tuple<2>));
        memcpy(process_data + prefix_sum_process_size[i], (&process_data_vector[i])->buffer, (&process_data_vector[i])->size);
        vector_buffer_free(&process_data_vector[i]);
    }
    free (process_data_vector);

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

#if 1

    for (u32 i3 = 0; i3 < g_t_actual_bucket_count; i3++)
    {
        u32 i = g_t_bucket_indices[i3];

        Relation1Map::iterator iy = delT[i].begin();
        for(; iy != delT[i].end(); iy++)
            delete (iy->second);
        delT[i].clear();
    }

#if DEBUG
    if (rank == 0)
        std::cout << "Global Comm End" << std::endl;

    if (rank == 0)
        std::cout << "Insert Start" << std::endl;
#endif

    double c2 = MPI_Wtime();
    double global_comm_time = (c2-c1);

    double i1 = MPI_Wtime();
    u64 tduplicates = 0;
    int tcount = 0;
    u64 count = 0;
    for (u32 i = 0; i < outer_hash_buffer_size; i = i + 2)
    {
        uint64_t bucket_id = hash_function(hash_buffer[i + 1]) % buckets;
        uint64_t sub_bucket_id = hash_function(hash_buffer[i]) % subbuckets_T[bucket_id];

        auto it = T[bucket_id].find(hash_buffer[i]);
        if ( it != T[bucket_id].end() ) {
            auto it2 = (it->second)->find(hash_buffer[i + 1]);
            if ( it2 != (it->second)->end() ) {
                tduplicates++;
            }
            else{
                (it->second)->insert(std::make_pair(hash_buffer[i + 1], 0));
                T[bucket_id].insert(std::make_pair(hash_buffer[i], it->second));
                tmap_sub_bucket[bucket_id][sub_bucket_id]++;

                tcount++;
                auto itx = delT[bucket_id].find(hash_buffer[i]);
                if ( itx != delT[bucket_id].end() ) {
                    auto it2x = (itx->second)->find(hash_buffer[i + 1]);
                    if ( it2x != (itx->second)->end() ) {
                        ;
                    }
                    else{
                        (itx->second)->insert(std::make_pair(hash_buffer[i + 1], 0));
                        delT[bucket_id].insert(std::make_pair(hash_buffer[i], itx->second));
                        count++;
                    }
                }
                else {
                    Relation0Map *k = new Relation0Map;
                    k->insert(std::make_pair(hash_buffer[i + 1], 0));
                    delT[bucket_id].insert(std::make_pair(hash_buffer[i],k));
                    count++;
                }
            }
        }
        else {
            Relation0Map *k = new Relation0Map;
            k->insert(std::make_pair(hash_buffer[i + 1], 0));
            T[bucket_id].insert(std::make_pair(hash_buffer[i],k));
            tmap_sub_bucket[bucket_id][sub_bucket_id]++;
            tcount++;

            auto itx = delT[bucket_id].find(hash_buffer[i]);
            if ( itx != delT[bucket_id].end() ) {
                auto it2x = (itx->second)->find(hash_buffer[i + 1]);
                if ( it2x != (itx->second)->end() ) {
                    ;
                }
                else{
                    (itx->second)->insert(std::make_pair(hash_buffer[i + 1], 0));
                    delT[bucket_id].insert(std::make_pair(hash_buffer[i], itx->second));
                    count++;
                }
            }
            else {
                Relation0Map *k = new Relation0Map;
                k->insert(std::make_pair(hash_buffer[i + 1], 0));
                delT[bucket_id].insert(std::make_pair(hash_buffer[i],k));
                count++;
            }
        }
    }
#endif
#endif

    delete[] hash_buffer;
    delete[] process_data;
    double i2 = MPI_Wtime();

#if DEBUG
    if (rank == 0)
        std::cout << "Insert End" << std::endl;
#endif

    double insert_time = (i2 - i1);

#if 1


    double v1 = MPI_Wtime();
    if (lc % 10 == 0)
    {
        int sum = 0;
        MPI_Allreduce(&count, &sum, 1, MPI_INT, MPI_BOR, comm);
        if(sum == 0)
            *lb = 1;
        else
            *lb = 0;
    }

    //
    *running_t_count = *running_t_count + tcount;
    double v2 = MPI_Wtime();


    *running_time = *running_time + (v2 - t1_s);

    if (rank == 0)
        std::cout << lc << " [" << v2-t1_s << "][" << c+jo+global_comm_time+insert_time+(v2-v1) << "][" << *running_time << "]"
                  << " N t " << tuple_count
                     //<< " non dup: " << non_deduplicate_tuple_count
                  << " Intra: " << c
                  << " Join: " << jo
                  << " C+J: " << c_j << " [" << c + jo << "]"
                  << " Comm: " << global_comm_time
                  << " Ins: " << insert_time
                  << " Ver: " << (v2 - v1)
                  << " Delta: " << tcount
                  << " Dup Delta: " << tduplicates
                  << " T : " << *running_t_count
                  << std::endl;

#endif
#endif
}


void load_balance_G(u32 buckets, u32*& gmap_bucket, u32**& gmap_sub_bucket, u32*& subbuckets_G, Relation1Map*& G, int**& gmap_sub_bucket_rank, int*& gmap_distinct_sub_bucket_rank_count, int**& gmap_distinct_sub_bucket_rank, int round_robin, float lb_factor)
{

#if DEBUG
    u64 G_tuple_count_before = 0;
    u64 G_tuple_count_before_global = 0;
    for (u32 bk = 0; bk < buckets; bk++)
    {
        for (u32 j = 0; j < subbuckets_G[bk]; j++)
            G_tuple_count_before = G_tuple_count_before + gmap_sub_bucket[bk][j];
    }
    MPI_Allreduce(&G_tuple_count_before, &G_tuple_count_before_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0)
        std::cout << "[G] Before Load balancing: " << G_tuple_count_before_global << std::endl;
#endif


    u32 max_sub_bucket_size = 0;
    u32 min_sub_bucket_size = INT_MAX;
    for (u32 i = 0; i < buckets; i++)
    {
        for (u32 j = 0; j < subbuckets_G[i]; j++)
        {
            if (gmap_sub_bucket[i][j] != 0)
            {
                if (gmap_sub_bucket[i][j] > max_sub_bucket_size)
                    max_sub_bucket_size = gmap_sub_bucket[i][j];

                if (gmap_sub_bucket[i][j] < min_sub_bucket_size)
                    min_sub_bucket_size = gmap_sub_bucket[i][j];
            }
        }
    }

    int global_max;
    int global_min;
    //u64 total_outer_hash_buffer_size = 0;
    MPI_Allreduce(&max_sub_bucket_size, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&min_sub_bucket_size, &global_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    u32 g_new_sub_bucket[buckets];
    memset(g_new_sub_bucket, 0, buckets * sizeof(u32));

    u32 global_g_new_sub_bucket[buckets];
    memset(global_g_new_sub_bucket, 0, buckets * sizeof(u32));
    for (u32 i = 0; i < buckets; i++)
        for (u32 j = 0; j < subbuckets_G[i]; j++)
            g_new_sub_bucket[i] = g_new_sub_bucket[i] + (gmap_sub_bucket[i][j] / (global_min*lb_factor));


    MPI_Allreduce(g_new_sub_bucket, global_g_new_sub_bucket, buckets, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    for (u32 i = 0; i < buckets; i++)
        if (global_g_new_sub_bucket[i] == 0)
            global_g_new_sub_bucket[i] = 1;

#if DEBUG
    if (rank == 0)
        std::cout << "[G Load Balance] New sub-bucket counts" << std::endl;
#endif

    int rcount = gmap_sub_bucket_rank[buckets-1][subbuckets_G[buckets-1] - 1] + 1;
    if (round_robin == 1)
    {
        for (u64 b = 0; b < buckets; b++)
        {
            if (subbuckets_G[b] < global_g_new_sub_bucket[b])
            {
                //gmap_bucket[b] = 0;
                std::unordered_set<int> distinct_g_ranks;

                int* temp = new int[subbuckets_G[b]];
                memcpy(temp, gmap_sub_bucket_rank[b], sizeof(int) * subbuckets_G[b]);

                delete[] gmap_sub_bucket_rank[b];
                gmap_sub_bucket_rank[b] = new int[global_g_new_sub_bucket[b]];

                memcpy(gmap_sub_bucket_rank[b], temp, sizeof(int) * subbuckets_G[b]);
                delete[] temp;

                for (u64 x = 0; x < global_g_new_sub_bucket[b]; x++)
                {
                    if (x >= subbuckets_G[b])
                        gmap_sub_bucket_rank[b][x] = rcount % nprocs;

                    if (gmap_sub_bucket_rank[b][x] == rank)
                        gmap_bucket[b] = 1;

                    distinct_g_ranks.insert(gmap_sub_bucket_rank[b][x]);
                    rcount++;
                }

                delete[] gmap_distinct_sub_bucket_rank[b];
                gmap_distinct_sub_bucket_rank[b] = new int[distinct_g_ranks.size()];
                u32 x = 0;
                for ( auto it = distinct_g_ranks.begin(); it != distinct_g_ranks.end(); ++it )
                {
                    gmap_distinct_sub_bucket_rank[b][x] = *it;
                    x++;
                }
                gmap_distinct_sub_bucket_rank_count[b] = x;
            }
        }

    }
    else
    {
        for (u64 b = 0; b < buckets; b++)
        {
            if (global_g_new_sub_bucket[b] != subbuckets_G[b])
            {
                //gmap_bucket[b] = 0;
                std::unordered_set<int> distinct_g_ranks;

                delete[] gmap_sub_bucket_rank[b];
                gmap_sub_bucket_rank[b] = new int[global_g_new_sub_bucket[b]];

                for (u64 x = 0; x < global_g_new_sub_bucket[b]; x++)
                {
                    gmap_sub_bucket_rank[b][x] = (hash_function((b << 32)^x))%nprocs;//rcount % nprocs;

                    if (gmap_sub_bucket_rank[b][x] == rank)
                        gmap_bucket[b] = 1;

                    distinct_g_ranks.insert(gmap_sub_bucket_rank[b][x]);
                }

                delete[] gmap_distinct_sub_bucket_rank[b];
                gmap_distinct_sub_bucket_rank[b] = new int[distinct_g_ranks.size()];
                u32 x = 0;
                for ( auto it = distinct_g_ranks.begin(); it != distinct_g_ranks.end(); ++it )
                {
                    gmap_distinct_sub_bucket_rank[b][x] = *it;
                    x++;
                }
                gmap_distinct_sub_bucket_rank_count[b] = x;
            }
        }
    }

#if DEBUG
    if (rank == 0)
        std::cout << "[G Load Balance] Rank allocation" << std::endl;
#endif

    /* process_size[j] stores the number of samples to be sent to process with rank j */
    int process_size[nprocs];
    memset(process_size, 0, nprocs * sizeof(int));

    /* vector[i] contains the data that needs to be sent to process i */
    vector_buffer* process_data_vector = (vector_buffer*)malloc(sizeof(vector_buffer) * nprocs);
    for (u32 i = 0; i < nprocs; ++i) {
        process_data_vector[i] = vector_buffer_create_empty();
    }

    int prefix_sum_process_size[nprocs];
    memset(prefix_sum_process_size, 0, nprocs * sizeof(int));

    int recv_process_size_buffer[nprocs];
    memset(recv_process_size_buffer, 0, nprocs * sizeof(int));


    int process_data_buffer_size = 0;
    u64 val[2] = {0,0};
    for (u32 i = 0; i < buckets; i++)
    {
        if (global_g_new_sub_bucket[i] != subbuckets_G[i])
        {
            delete[] gmap_sub_bucket[i];
            gmap_sub_bucket[i] = new u32[global_g_new_sub_bucket[i]];
            memset(gmap_sub_bucket[i], 0, sizeof(u32) * global_g_new_sub_bucket[i]);


            for (auto it = G[i].begin(); it != G[i].end(); it++)
            {
                Relation0Map* it2 = it->second;
                for (auto dit2 = it2->begin(); dit2 != it2->end(); dit2++)
                {
                    uint64_t bucket_id = hash_function(it->first) % buckets;
                    uint64_t sub_bucket_id = hash_function(dit2->first) % global_g_new_sub_bucket[bucket_id];

                    //uint64_t index = (hash_function((bucket_id << 32)^sub_bucket_id))%nprocs;
                    int index = gmap_sub_bucket_rank[bucket_id][sub_bucket_id];
                    process_size[index] = process_size[index] + COL_COUNT;

                    val[0] = it->first;
                    val[1] = dit2->first;
                    vector_buffer_append(&process_data_vector[index],(unsigned char *) val, sizeof(u64)*2);
                }
            }

            Relation1Map::iterator iy2 = G[i].begin();
            for(; iy2 != G[i].end(); iy2++)
                delete (iy2->second);

            G[i].clear();
        }
    }

#if DEBUG
    if (rank == 0)
        std::cout << "[G Load Balance] Collecting data to send" << std::endl;
#endif

    for (u32 ix = 1; ix < nprocs; ix++)
        prefix_sum_process_size[ix] = prefix_sum_process_size[ix - 1] + process_size[ix - 1];

    process_data_buffer_size = prefix_sum_process_size[nprocs - 1] + process_size[nprocs - 1];

    u64* process_data = new u64[process_data_buffer_size];

    for (u32 ix = 0; ix < nprocs; ix++)
    {
        memcpy(process_data + prefix_sum_process_size[ix], (&process_data_vector[ix])->buffer, (&process_data_vector[ix])->size);
        vector_buffer_free(&process_data_vector[ix]);
    }
    free (process_data_vector);


    MPI_Alltoall(process_size, 1, MPI_INT, recv_process_size_buffer, 1, MPI_INT, MPI_COMM_WORLD);

#if 1
    int prefix_sum_recv_process_size_buffer[nprocs];
    memset(prefix_sum_recv_process_size_buffer, 0, nprocs * sizeof(int));
    for (u32 ix = 1; ix < nprocs; ix++)
        prefix_sum_recv_process_size_buffer[ix] = prefix_sum_recv_process_size_buffer[ix - 1] + recv_process_size_buffer[ix - 1];


    int outer_hash_buffer_size = 0;
    for(u32 ix = 0; ix < nprocs; ix++)
        outer_hash_buffer_size = outer_hash_buffer_size + recv_process_size_buffer[ix];


    u64* outer_hash_data = new u64[outer_hash_buffer_size];
    MPI_Alltoallv(process_data, process_size, prefix_sum_process_size, MPI_UNSIGNED_LONG_LONG, outer_hash_data, recv_process_size_buffer, prefix_sum_recv_process_size_buffer, MPI_UNSIGNED_LONG_LONG, comm);


    delete[] process_data;

#if DEBUG
    if (rank == 0)
        std::cout << "[G Load Balance] All to All" << std::endl;
#endif
    //subbuckets_G[i] = global_g_new_sub_bucket[i];
    for (int in = 0; in < outer_hash_buffer_size; in = in + 2)
    {
        uint64_t bucket_id = hash_function((outer_hash_data)[in]) % buckets;
        uint64_t sub_bucket_id = hash_function((outer_hash_data)[in+1]) % global_g_new_sub_bucket[bucket_id];
        gmap_bucket[bucket_id] = 1;


        //if (rank == 1)
        //std::cout << "Val " << outer_hash_data[in] << " b " << bucket_id << " sbi " << sub_bucket_id << std::endl;

        auto it = G[bucket_id].find((outer_hash_data)[in]);
        if( it != G[bucket_id].end() ) {
            auto it2 = (it->second)->find((outer_hash_data)[in+1]);
            if( it2 != (it->second)->end() ) {
                ;
            }
            else{
                (it->second)->insert(std::make_pair((outer_hash_data)[in + 1], 0));
                G[bucket_id][(outer_hash_data)[in]] = it->second;
                gmap_sub_bucket[bucket_id][sub_bucket_id]++;
            }
        }
        else {
            Relation0Map *k = new Relation0Map;
            k->insert(std::make_pair((outer_hash_data)[in + 1], 0));
            G[bucket_id].insert(std::make_pair((outer_hash_data)[in],k));
            gmap_sub_bucket[bucket_id][sub_bucket_id]++;
        }
    }

    delete[] outer_hash_data;
#endif

    if (rank == 0)
    {
        for (u32 i = 0; i < buckets; i++)
            std::cout << i << "[G] Bucket count changed from " << subbuckets_G[i] << " to " << global_g_new_sub_bucket[i] << std::endl;
    }
    memcpy(subbuckets_G, global_g_new_sub_bucket, sizeof(u32)*buckets);

#if DEBUG
    if (rank == 0)
        std::cout << "[G Load Balance] Local inserts" << std::endl;
#endif

#if DEBUG
    u64 G_tuple_count_after = 0;
    u64 G_tuple_count_after_global = 0;
    for (u32 bk = 0; bk < buckets; bk++)
    {
        for (u32 j = 0; j < subbuckets_G[bk]; j++)
            G_tuple_count_after = G_tuple_count_after + gmap_sub_bucket[bk][j];
    }
    MPI_Allreduce(&G_tuple_count_after, &G_tuple_count_after_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    assert(G_tuple_count_before_global == G_tuple_count_after_global);
    if (rank == 0)
        std::cout << "[G] After Load balanging " << G_tuple_count_before_global << " " << G_tuple_count_after_global << std::endl;

    std::cout << "Balanced [G] Rank: " << rank << " Hashed count: " << G_tuple_count_after << " Ratio: " << (float)((float) (G_tuple_count_after)/G_tuple_count_after_global)*100 << std::endl;
#endif

}



void load_balance_T(u32 buckets, u32*& tmap_bucket, u32**& tmap_sub_bucket, u32*& subbuckets_T, Relation1Map*& T, Relation1Map*& delta, int**& tmap_sub_bucket_rank, int*& tmap_distinct_sub_bucket_rank_count, int**& tmap_distinct_sub_bucket_rank, int round_robin, float lb_factor)
{
#if DEBUG
    u64 T_tuple_count_before = 0;
    u64 T_tuple_count_before_global = 0;

    for (u32 bk = 0; bk < buckets; bk++)
    {
        for (u32 j = 0; j < subbuckets_T[bk]; j++)
            T_tuple_count_before = T_tuple_count_before + tmap_sub_bucket[bk][j];
    }
    MPI_Allreduce(&T_tuple_count_before, &T_tuple_count_before_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0)
        std::cout << "[T] Before Load balancing " << T_tuple_count_before_global << std::endl;
#endif


    u32 max_sub_bucket_size = 0;
    u32 min_sub_bucket_size = INT_MAX;


    for (u32 i = 0; i < buckets; i++)
    {
        for (u32 j = 0; j < subbuckets_T[i]; j++)
        {
            if (tmap_sub_bucket[i][j] != 0)
            {
                if (tmap_sub_bucket[i][j] > max_sub_bucket_size)
                    max_sub_bucket_size = tmap_sub_bucket[i][j];

                if (tmap_sub_bucket[i][j] < min_sub_bucket_size)
                    min_sub_bucket_size = tmap_sub_bucket[i][j];
            }
        }
    }

    int global_max;
    int global_min;
    MPI_Allreduce(&max_sub_bucket_size, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&min_sub_bucket_size, &global_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    u32 t_new_sub_bucket[buckets];
    memset(t_new_sub_bucket, 0, buckets * sizeof(u32));

    u32 global_t_new_sub_bucket[buckets];
    memset(global_t_new_sub_bucket, 0, buckets * sizeof(u32));
    for (u32 i = 0; i < buckets; i++)
        for (u32 j = 0; j < subbuckets_T[i]; j++)
            t_new_sub_bucket[i] = t_new_sub_bucket[i] + (tmap_sub_bucket[i][j] / (global_min*lb_factor) );

    MPI_Allreduce(t_new_sub_bucket, global_t_new_sub_bucket, buckets, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    u32 ct = 0;
    for (u32 i = 0; i < buckets; i++)
    {
        if (global_t_new_sub_bucket[i] == 0)
            global_t_new_sub_bucket[i] = 1;
        if (global_t_new_sub_bucket[i] > nprocs/4)
            global_t_new_sub_bucket[i] = nprocs/4;

        if (global_t_new_sub_bucket[i] == subbuckets_T[i])
            ct++;
    }
    if (ct == buckets)
        return;


#if DEBUG
    if (rank == 0)
        std::cout << "[T Load Balance] New sub-bucket counts" << std::endl;
#endif


    int rcount = tmap_sub_bucket_rank[buckets-1][subbuckets_T[buckets-1] - 1] + 1;
    //std::cout << "rcount is " << rcount << std::endl;
    if (round_robin == 1)
    {
        for (u64 b = 0; b < buckets; b++)
        {
            //assert(subbuckets_T[b] < global_t_new_sub_bucket[b]);
            if (subbuckets_T[b] < global_t_new_sub_bucket[b])
            {
                //tmap_bucket[b] = 0;

                std::unordered_set<int> distinct_t_ranks;
                int* temp = new int[subbuckets_T[b]];
                memcpy(temp, tmap_sub_bucket_rank[b], sizeof(int) * subbuckets_T[b]);

                delete[] tmap_sub_bucket_rank[b];
                tmap_sub_bucket_rank[b] = new int[global_t_new_sub_bucket[b]];

                memcpy(tmap_sub_bucket_rank[b], temp, sizeof(int) * subbuckets_T[b]);
                delete[] temp;

                for (u64 x = 0; x < global_t_new_sub_bucket[b]; x++)
                {
                    if (x >= subbuckets_T[b])
                        tmap_sub_bucket_rank[b][x] = rcount % nprocs;

                    if (tmap_sub_bucket_rank[b][x] == rank)
                        tmap_bucket[b] = 1;

                    distinct_t_ranks.insert(tmap_sub_bucket_rank[b][x]);
                    rcount++;
                }

                delete[] tmap_distinct_sub_bucket_rank[b];
                tmap_distinct_sub_bucket_rank[b] = new int[distinct_t_ranks.size()];
                u32 x = 0;
                for ( auto it = distinct_t_ranks.begin(); it != distinct_t_ranks.end(); ++it )
                {
                    tmap_distinct_sub_bucket_rank[b][x] = *it;
                    x++;
                }
                tmap_distinct_sub_bucket_rank_count[b] = x;
            }
            else if (subbuckets_T[b] > global_t_new_sub_bucket[b])
            {
                std::unordered_set<int> distinct_t_ranks;
                int* temp = new int[subbuckets_T[b]];
                memcpy(temp, tmap_sub_bucket_rank[b], sizeof(int) * subbuckets_T[b]);

                delete[] tmap_sub_bucket_rank[b];
                tmap_sub_bucket_rank[b] = new int[global_t_new_sub_bucket[b]];

                memcpy(tmap_sub_bucket_rank[b], temp, sizeof(int) * global_t_new_sub_bucket[b]);
                delete[] temp;

                for (u64 x = 0; x < global_t_new_sub_bucket[b]; x++)
                {
                    //if (x >= subbuckets_T[b])
                    //    tmap_sub_bucket_rank[b][x] = rcount % nprocs;

                    if (tmap_sub_bucket_rank[b][x] == rank)
                        tmap_bucket[b] = 1;

                    distinct_t_ranks.insert(tmap_sub_bucket_rank[b][x]);
                    rcount++;
                }

                delete[] tmap_distinct_sub_bucket_rank[b];
                tmap_distinct_sub_bucket_rank[b] = new int[distinct_t_ranks.size()];
                u32 x = 0;
                for ( auto it = distinct_t_ranks.begin(); it != distinct_t_ranks.end(); ++it )
                {
                    tmap_distinct_sub_bucket_rank[b][x] = *it;
                    x++;
                }
                tmap_distinct_sub_bucket_rank_count[b] = x;
            }
        }
    }
    else
    {
        for (u64 b = 0; b < buckets; b++)
        {
            if (subbuckets_T[b] != global_t_new_sub_bucket[b])
            {
                tmap_bucket[b] = 0;

                std::unordered_set<int> distinct_t_ranks;

                delete[] tmap_sub_bucket_rank[b];
                tmap_sub_bucket_rank[b] = new int[global_t_new_sub_bucket[b]];


                //delete[] tmap_sub_bucket_rank[b];
                //tmap_sub_bucket_rank[b] = new int[global_t_new_sub_bucket[b]];

                for (u64 x = 0; x < global_t_new_sub_bucket[b]; x++)
                {
                    //if (x >= subbuckets_T[b])
                    tmap_sub_bucket_rank[b][x] = (hash_function((b << 32)^x))%nprocs;//rcount % nprocs;

                    if (tmap_sub_bucket_rank[b][x] == rank)
                        tmap_bucket[b] = 1;

                    distinct_t_ranks.insert(tmap_sub_bucket_rank[b][x]);
                    rcount++;
                }

                delete[] tmap_distinct_sub_bucket_rank[b];
                tmap_distinct_sub_bucket_rank[b] = new int[distinct_t_ranks.size()];
                u32 x = 0;
                for ( auto it = distinct_t_ranks.begin(); it != distinct_t_ranks.end(); ++it )
                {
                    tmap_distinct_sub_bucket_rank[b][x] = *it;
                    x++;
                }
                tmap_distinct_sub_bucket_rank_count[b] = x;
            }
        }
    }

#if DEBUG
    if (rank == 0)
        std::cout << "[T Load Balance] Rank allocation" << std::endl;
#endif

#if 1

    /* process_size[j] stores the number of samples to be sent to process with rank j */
    int process_size[nprocs];
    memset(process_size, 0, nprocs * sizeof(int));

    /* vector[i] contains the data that needs to be sent to process i */
    vector_buffer* process_data_vector = (vector_buffer*)malloc(sizeof(vector_buffer) * nprocs);
    for (u32 i = 0; i < nprocs; ++i) {
        process_data_vector[i] = vector_buffer_create_empty();
    }

    int prefix_sum_process_size[nprocs];
    memset(prefix_sum_process_size, 0, nprocs * sizeof(int));

    int recv_process_size_buffer[nprocs];
    memset(recv_process_size_buffer, 0, nprocs * sizeof(int));

    u64 val[2] = {0, 0};
    for (u32 bk = 0; bk < buckets; bk++)
    {
        if (subbuckets_T[bk] != global_t_new_sub_bucket[bk])
        {
            delete[] tmap_sub_bucket[bk];
            tmap_sub_bucket[bk] = new u32[global_t_new_sub_bucket[bk]];
            memset(tmap_sub_bucket[bk], 0, sizeof(u32) * global_t_new_sub_bucket[bk]);

            for (auto it = T[bk].begin(); it != T[bk].end(); it++)
            {
                Relation0Map* it2 = it->second;
                for (auto dit2 = it2->begin(); dit2 != it2->end(); dit2++)
                {
                    //test
                    uint64_t bucket_id = hash_function(dit2->first) % buckets;
                    uint64_t sub_bucket_id = hash_function(it->first) % global_t_new_sub_bucket[bucket_id];

                    //uint64_t index = (hash_function((bucket_id << 32)^sub_bucket_id))%nprocs;
                    int index = tmap_sub_bucket_rank[bucket_id][sub_bucket_id];
                    process_size[index] = process_size[index] + COL_COUNT;

                    val[0] = it->first;
                    val[1] = dit2->first;
                    vector_buffer_append(&process_data_vector[index],(unsigned char *) val, sizeof(u64)*2);
                }
            }
            Relation1Map::iterator iy2 = T[bk].begin();
            for(; iy2 != T[bk].end(); iy2++)
                delete (iy2->second);

            T[bk].clear();
        }
    }

#if DEBUG
    if (rank == 0)
        std::cout << "[T Load Balance] Collecting data to send" << std::endl;
#endif

    for (u32 ix = 1; ix < nprocs; ix++)
        prefix_sum_process_size[ix] = prefix_sum_process_size[ix - 1] + process_size[ix - 1];

    int process_data_buffer_size = prefix_sum_process_size[nprocs - 1] + process_size[nprocs - 1];

    //u64 process_data[process_data_buffer_size];
    u64* process_data = new u64[process_data_buffer_size];

    for (u32 ix = 0; ix < nprocs; ix++)
    {
        memcpy(process_data + prefix_sum_process_size[ix], (&process_data_vector[ix])->buffer, (&process_data_vector[ix])->size);
        vector_buffer_free(&process_data_vector[ix]);
    }
    free(process_data_vector);


    /* This step prepares for actual data transfer */
    /* Every process sends to every other process the amount of data it is going to send */


    MPI_Alltoall(process_size, 1, MPI_INT, recv_process_size_buffer, 1, MPI_INT, MPI_COMM_WORLD);


    int prefix_sum_recv_process_size_buffer[nprocs];
    memset(prefix_sum_recv_process_size_buffer, 0, nprocs * sizeof(int));
    for (u32 ix = 1; ix < nprocs; ix++)
        prefix_sum_recv_process_size_buffer[ix] = prefix_sum_recv_process_size_buffer[ix - 1] + recv_process_size_buffer[ix - 1];


    /* Sending data to all processes */
    /* What is the buffer size to allocate */
    int outer_hash_buffer_size = 0;
    for(u32 ix = 0; ix < nprocs; ix++)
        outer_hash_buffer_size = outer_hash_buffer_size + recv_process_size_buffer[ix];


    u64* outer_hash_data = new u64[outer_hash_buffer_size];
    MPI_Alltoallv(process_data, process_size, prefix_sum_process_size, MPI_UNSIGNED_LONG_LONG, outer_hash_data, recv_process_size_buffer, prefix_sum_recv_process_size_buffer, MPI_UNSIGNED_LONG_LONG, comm);

    delete[] process_data;

#if DEBUG
    if (rank == 0)
        std::cout << "[T Load Balance] All to All" << std::endl;
#endif

#if 1
    //std::cout << "outer_hash_buffer_size " << outer_hash_buffer_size << std::endl;
    //subbuckets_G[i] = global_t_new_sub_bucket[i];
    for (int in = 0; in < outer_hash_buffer_size; in = in + 2)
    {
        uint64_t bucket_id = hash_function((outer_hash_data)[in + 1]) % buckets;
        uint64_t sub_bucket_id = hash_function((outer_hash_data)[in]) % global_t_new_sub_bucket[bucket_id];
        tmap_bucket[bucket_id] = 1;


        auto it = T[bucket_id].find((outer_hash_data)[in]);
        if( it != T[bucket_id].end() ) {
            auto it2 = (it->second)->find((outer_hash_data)[in+1]);
            if( it2 != (it->second)->end() ) {
                ;
            }
            else{
                (it->second)->insert(std::make_pair((outer_hash_data)[in + 1], 0));
                T[bucket_id][(outer_hash_data)[in]] = it->second;
                tmap_sub_bucket[bucket_id][sub_bucket_id]++;
            }
        }
        else {
            Relation0Map *k = new Relation0Map;
            k->insert(std::make_pair((outer_hash_data)[in + 1], 0));
            T[bucket_id].insert(std::make_pair((outer_hash_data)[in],k));
            tmap_sub_bucket[bucket_id][sub_bucket_id]++;
        }
    }

    delete[] outer_hash_data;

#if DEBUG
    if (rank == 0)
        std::cout << "[T Load Balance] Local inserts" << std::endl;
#endif

#if DEBUG
    u64 T_tuple_count_after = 0;
    u64 T_tuple_count_after_global = 0;
    for (u32 bk = 0; bk < buckets; bk++)
    {
        for (u32 j = 0; j < global_t_new_sub_bucket[bk]; j++)
            T_tuple_count_after = T_tuple_count_after + tmap_sub_bucket[bk][j];
    }
    MPI_Allreduce(&T_tuple_count_after, &T_tuple_count_after_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0)
        std::cout << "[T] After Load balancing " << T_tuple_count_before_global << " " << T_tuple_count_after_global << std::endl;
    assert(T_tuple_count_before_global == T_tuple_count_after_global);
#endif


#endif


#endif

#if 1
#if DEBUG
    if (rank == 0)
        std::cout << "[dT] Before Load balancing " << std::endl;
#endif

    /* process_size[j] stores the number of samples to be sent to process with rank j */
    int process_size_dt[nprocs];
    memset(process_size_dt, 0, nprocs * sizeof(int));

    /* vector[i] contains the data that needs to be sent to process i */
    vector_buffer* process_data_vector_dt1 = (vector_buffer*)malloc(sizeof(vector_buffer) * nprocs);
    for (u32 i = 0; i < nprocs; ++i) {
        process_data_vector_dt1[i] = vector_buffer_create_empty();
    }

    int prefix_sum_process_size_dt[nprocs];
    memset(prefix_sum_process_size_dt, 0, nprocs * sizeof(int));

    int recv_process_size_buffer_dt[nprocs];
    memset(recv_process_size_buffer_dt, 0, nprocs * sizeof(int));

    for (u32 bk = 0; bk < buckets; bk++)
    {
        if (subbuckets_T[bk] != global_t_new_sub_bucket[bk])
        {
            for (auto it = delta[bk].begin(); it != delta[bk].end(); it++)
            {
                Relation0Map* it2 = it->second;
                for (auto dit2 = it2->begin(); dit2 != it2->end(); dit2++)
                {
                    uint64_t bucket_id = hash_function(dit2->first) % buckets;
                    uint64_t sub_bucket_id = hash_function(it->first) % global_t_new_sub_bucket[bucket_id];

                    int index = tmap_sub_bucket_rank[bucket_id][sub_bucket_id];
                    //uint64_t index = (hash_function((bucket_id << 32)^sub_bucket_id))%nprocs;
                    process_size_dt[index] = process_size_dt[index] + COL_COUNT;

                    val[0] = it->first;
                    val[1] = dit2->first;
                    vector_buffer_append(&process_data_vector_dt1[index],(unsigned char *) val, sizeof(u64)*2);
                }
            }

            Relation1Map::iterator iy2 = delta[bk].begin();
            for(; iy2 != delta[bk].end(); iy2++)
                delete (iy2->second);

            delta[bk].clear();
        }
    }
#if DEBUG
    if (rank == 0)
        std::cout << "[dT Load Balance] Collecting data to send" << std::endl;
#endif

    for (u32 ix = 1; ix < nprocs; ix++)
        prefix_sum_process_size_dt[ix] = prefix_sum_process_size_dt[ix - 1] + process_size_dt[ix - 1];

    int process_data_buffer_size_dt = prefix_sum_process_size_dt[nprocs - 1] + process_size_dt[nprocs - 1];

    //u64 process_data[process_data_buffer_size];
    u64* process_data_dt = new u64[process_data_buffer_size_dt];

    for (u32 ix = 0; ix < nprocs; ix++)
    {
        memcpy(process_data_dt + prefix_sum_process_size_dt[ix], (&process_data_vector_dt1[ix])->buffer, (&process_data_vector_dt1[ix])->size);
        vector_buffer_free(&process_data_vector_dt1[ix]);
    }
    free (process_data_vector_dt1);


    /* This step prepares for actual data transfer */
    /* Every process sends to every other process the amount of data it is going to send */


    MPI_Alltoall(process_size_dt, 1, MPI_INT, recv_process_size_buffer_dt, 1, MPI_INT, MPI_COMM_WORLD);

#if 1
    int prefix_sum_recv_process_size_buffer_dt[nprocs];
    memset(prefix_sum_recv_process_size_buffer_dt, 0, nprocs * sizeof(int));
    for (u32 ix = 1; ix < nprocs; ix++)
        prefix_sum_recv_process_size_buffer_dt[ix] = prefix_sum_recv_process_size_buffer_dt[ix - 1] + recv_process_size_buffer_dt[ix - 1];


    /* Sending data to all processes */
    /* What is the buffer size to allocate */
    int outer_hash_buffer_size_dt = 0;
    for (u32 ix = 0; ix < nprocs; ix++)
        outer_hash_buffer_size_dt = outer_hash_buffer_size_dt + recv_process_size_buffer_dt[ix];


    u64* outer_hash_data_dt = new u64[outer_hash_buffer_size_dt];
    MPI_Alltoallv(process_data_dt, process_size_dt, prefix_sum_process_size_dt, MPI_UNSIGNED_LONG_LONG, outer_hash_data_dt, recv_process_size_buffer_dt, prefix_sum_recv_process_size_buffer_dt, MPI_UNSIGNED_LONG_LONG, comm);


    delete[] process_data_dt;

#if DEBUG
    if (rank == 0)
        std::cout << "[dT Load Balance] All to All" << std::endl;
#endif


    //std::cout << "outer_hash_buffer_size " << outer_hash_buffer_size << std::endl;
    //subbuckets_G[i] = global_t_new_sub_bucket[i];
    for (int in = 0; in < outer_hash_buffer_size_dt; in = in + 2)
    {
        uint64_t bucket_id = hash_function((outer_hash_data_dt)[in+1]) % buckets;

        auto it = delta[bucket_id].find((outer_hash_data_dt)[in]);
        if( it != delta[bucket_id].end() ) {
            auto it2 = (it->second)->find((outer_hash_data_dt)[in+1]);
            if( it2 != (it->second)->end() ) {
                ;
            }
            else{
                (it->second)->insert(std::make_pair((outer_hash_data_dt)[in + 1], 0));
                delta[bucket_id].insert(std::make_pair((outer_hash_data_dt)[in], it->second));
            }
        }
        else {
            Relation0Map *k = new Relation0Map;
            k->insert(std::make_pair((outer_hash_data_dt)[in + 1], 0));
            delta[bucket_id].insert(std::make_pair((outer_hash_data_dt)[in],k));
        }
    }
#endif
    delete[] outer_hash_data_dt;

#if DEBUG
    if (rank == 0)
        std::cout << "[dT Load Balance] Local inserts" << std::endl;
#endif

    if (rank == 0)
        for (u32 i = 0; i < buckets; i++)
            std::cout << "[T] " << i <<" Bucket count changed from " << subbuckets_T[i] << " to " << global_t_new_sub_bucket[i] << std::endl;

    memcpy(subbuckets_T, global_t_new_sub_bucket, sizeof(u32)*buckets);

    //std::cout << "Balanced [T] Rank: " << rank << " Hashed count: " << outer_hash_buffer_size << " Ratio: " << (float)((float) outer_hash_buffer_size/global_row_count)*100 << std::endl;
#if DEBUG
    if (rank == 0)
        std::cout << "[dT] After Load balancing " << std::endl;


    assert(T_tuple_count_before_global == T_tuple_count_after_global);

    std::cout << "Balanced [T] Rank: " << rank << " Hashed count: " << T_tuple_count_after << " Ratio: " << (float)((float) T_tuple_count_after/T_tuple_count_after_global)*100 << std::endl;
#endif
#endif

}



void dump_relation(Relation1Map*** X, int lc, int buckets, u32* subbuckets, u32* map)
{

    for (int i = 0; i < buckets; i++)
    {

        if (map[i] == 1)
        {
            std::ofstream myfile;
            char filename[1024];
            sprintf(filename, "filename_%d_%d_%d", lc, rank, i);
            myfile.open (filename);

            for (u32 j = 0; j < subbuckets[i]; j++)
            {
                for ( auto local_it = X[i][j]->begin(); local_it!= X[i][j]->end(); ++local_it )
                {
                    Relation0Map* k = local_it->second;
                    for (auto it2 = k->begin(); it2 != k->end(); it2++)
                        myfile << local_it->first << " " << it2->first << std::endl;
                }
            }

            myfile.close();
        }
    }
}



int main(int argc, char **argv)
{
    // Initializing MPI
    MPI_Init(&argc, &argv);
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc != 7)
    {
        if (rank == 0)
        {
            std::cout << "Usage:" << "TC <bucket count factor> <initial sub bucket count> <comm type> <balance type> <load balancing factor>"   << std::endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }

#if 1
    comm = MPI_COMM_WORLD;
    double ior_start = MPI_Wtime();

    nprocs = (u32)size;

    u32 buckets;
    u32 *subbuckets_G;
    u32 *subbuckets_T;

    // this is the initialization step, we can change this as well
    buckets = nprocs * atoi(argv[2]);

    // the relations can have different number of sub buckets
    subbuckets_G = new u32[buckets];
    memset(subbuckets_G, 0, buckets*sizeof(u32));

    subbuckets_T = new u32[buckets];
    memset(subbuckets_T, 0, buckets*sizeof(u32));

    for (u64 b = 0; b < buckets; b++)
    {
        subbuckets_G[b] = atoi(argv[3]);
        subbuckets_T[b] = atoi(argv[3]);
    }

    int** gmap_sub_bucket_rank = new int*[buckets];
    for (u64 b = 0; b < buckets; b++)
        gmap_sub_bucket_rank[b] = new int[subbuckets_G[b]];

    int** tmap_sub_bucket_rank = new int*[buckets];
    for (u64 b = 0; b < buckets; b++)
        tmap_sub_bucket_rank[b] = new int[subbuckets_T[b]];


    int** gmap_distinct_sub_bucket_rank = new int*[buckets];
    int* gmap_distinct_sub_bucket_rank_count = new int[buckets];


    u32 *gmap_bucket = new u32[buckets];
    memset(gmap_bucket, 0, sizeof(u32) * buckets);

    int** tmap_distinct_sub_bucket_rank = new int*[buckets];
    int* tmap_distinct_sub_bucket_rank_count = new int[buckets];

    u32 *tmap_bucket = new u32[buckets];
    memset(tmap_bucket, 0, sizeof(u32) * buckets);


    int rcount = 0;

    // hashing based distribution
    if (atoi(argv[4]) == 0)
    {
        for (u64 b = 0; b < buckets; b++)
        {
            std::unordered_set<int> distinct_g_ranks;
            for (u64 x = 0; x < subbuckets_G[b]; x++)
            {
                gmap_sub_bucket_rank[b][x] = (hash_function((b << 32)^x))%nprocs;//rcount % nprocs;

                if (gmap_sub_bucket_rank[b][x] == rank)
                    gmap_bucket[b] = 1;

                distinct_g_ranks.insert(gmap_sub_bucket_rank[b][x]);
                rcount++;
            }
            gmap_distinct_sub_bucket_rank[b] = new int[distinct_g_ranks.size()];
            u32 x = 0;
            for ( auto it = distinct_g_ranks.begin(); it != distinct_g_ranks.end(); ++it )
            {
                gmap_distinct_sub_bucket_rank[b][x] = *it;
                x++;
            }
            gmap_distinct_sub_bucket_rank_count[b] = x;
        }

        rcount = 0;
        for (u64 b = 0; b < buckets; b++)
        {
            std::unordered_set<int> distinct_t_ranks;
            for (u64 x = 0; x < subbuckets_T[b]; x++)
            {
                tmap_sub_bucket_rank[b][x] = (hash_function((b << 32)^x))%nprocs;//rcount % nprocs;

                if (tmap_sub_bucket_rank[b][x] == rank)
                    tmap_bucket[b] = 1;

                distinct_t_ranks.insert(tmap_sub_bucket_rank[b][x]);
                rcount++;
            }
            tmap_distinct_sub_bucket_rank[b] = new int[distinct_t_ranks.size()];
            u32 x = 0;
            for ( auto it = distinct_t_ranks.begin(); it != distinct_t_ranks.end(); ++it )
            {
                tmap_distinct_sub_bucket_rank[b][x] = *it;
                x++;
            }
            tmap_distinct_sub_bucket_rank_count[b] = x;
        }
    }
    else
    {
        rcount = 0;
        for (u64 b = 0; b < buckets; b++)
        {
            std::unordered_set<int> distinct_g_ranks;
            for (u64 x = 0; x < subbuckets_G[b]; x++)
            {
                gmap_sub_bucket_rank[b][x] = rcount % nprocs;

                if (gmap_sub_bucket_rank[b][x] == rank)
                    gmap_bucket[b] = 1;

                distinct_g_ranks.insert(gmap_sub_bucket_rank[b][x]);
                rcount++;
            }
            gmap_distinct_sub_bucket_rank[b] = new int[distinct_g_ranks.size()];
            u32 x = 0;
            for ( auto it = distinct_g_ranks.begin(); it != distinct_g_ranks.end(); ++it )
            {
                gmap_distinct_sub_bucket_rank[b][x] = *it;
                x++;
            }
            gmap_distinct_sub_bucket_rank_count[b] = x;
        }

        rcount = 0;
        for (u64 b = 0; b < buckets; b++)
        {
            std::unordered_set<int> distinct_t_ranks;
            for (u64 x = 0; x < subbuckets_T[b]; x++)
            {
                tmap_sub_bucket_rank[b][x] = rcount % nprocs;

                if (tmap_sub_bucket_rank[b][x] == rank)
                    tmap_bucket[b] = 1;

                distinct_t_ranks.insert(tmap_sub_bucket_rank[b][x]);
                rcount++;
            }
            tmap_distinct_sub_bucket_rank[b] = new int[distinct_t_ranks.size()];
            u32 x = 0;
            for ( auto it = distinct_t_ranks.begin(); it != distinct_t_ranks.end(); ++it )
            {
                tmap_distinct_sub_bucket_rank[b][x] = *it;
                x++;
            }
            tmap_distinct_sub_bucket_rank_count[b] = x;
        }
    }


    u32 entry_count;
    u64 *input_buffer = NULL;
    parallel_read_input_relation_from_file_to_local_buffer(argv[1], &input_buffer, &entry_count);

    double ior_end = MPI_Wtime();

    double hash_start = MPI_Wtime();

    u32 G_hash_entry_count;
    u64 *G_hashed_data = NULL;
    buffer_data_to_hash_buffer(entry_count, input_buffer, 0, &G_hashed_data, &G_hash_entry_count, MPI_COMM_WORLD, buckets, subbuckets_G, subbuckets_T, gmap_sub_bucket_rank);

    std::cout << "[G] Rank: " << rank << " Hashed count: " << G_hash_entry_count << " Ratio: " << (float)((float) G_hash_entry_count/global_row_count)*100 << std::endl;


    u32 T_hash_entry_count;
    u64 *T_hashed_data = NULL;
    buffer_data_to_hash_buffer(entry_count, input_buffer,  1, &T_hashed_data, &T_hash_entry_count, MPI_COMM_WORLD, buckets, subbuckets_G, subbuckets_T, tmap_sub_bucket_rank);

    std::cout << "[T] Rank: " << rank << " Hashed count: " << T_hash_entry_count << " Ratio: " << (float)((float) T_hash_entry_count/global_row_count)*100 << std::endl;

    delete[] input_buffer;
    double hash_end = MPI_Wtime();

    double relation_start = MPI_Wtime();

    // this buffer is used as a bitmap

    u32 **gmap_sub_bucket = new u32*[buckets];
    memset(gmap_sub_bucket, 0, sizeof(u32*) * buckets);

    for (u32 i = 0; i < buckets; i++)
    {
        gmap_sub_bucket[i] = new u32[subbuckets_G[i]];
        memset(gmap_sub_bucket[i], 0, sizeof(u32) * subbuckets_G[i]);
    }

    Relation1Map *G = new Relation1Map[buckets];


    for (u32 i = 0; i < COL_COUNT * G_hash_entry_count; i = i + 2)
    {
        uint64_t bucket_id = hash_function(G_hashed_data[i]) % buckets;
        uint64_t sub_bucket_id = hash_function(G_hashed_data[i+1]) % subbuckets_G[bucket_id];
        gmap_bucket[bucket_id] = 1;


        auto it = G[bucket_id].find(G_hashed_data[i]);
        if( it != G[bucket_id].end() ) {
            auto it2 = (it->second)->find(G_hashed_data[i + 1]);
            if( it2 != (it->second)->end() ) {
                ;
            }
            else{
                (it->second)->insert(std::make_pair(G_hashed_data[i + 1], 0));
                G[bucket_id][G_hashed_data[i]] = it->second;
                gmap_sub_bucket[bucket_id][sub_bucket_id]++;
            }
        }
        else {
            Relation0Map *k = new Relation0Map;
            k->insert(std::make_pair(G_hashed_data[i + 1], 0));
            G[bucket_id].insert(std::make_pair(G_hashed_data[i],k));
            gmap_sub_bucket[bucket_id][sub_bucket_id]++;
        }
    }


#if 0
    for (u32 i = 0; i < buckets; i++)
    {
        for (u32 j = 0; j < subbuckets_G[i]; j++)
        {
            int counter = 0;
            for (auto it = G[i][j].begin(); it != G[i][j].end(); it++)
            {
                Relation0Map* it2 = it->second;
                for (auto dit2 = it2->begin(); dit2 != it2->end(); dit2++)
                {
                    //std::cout << "Rank " << rank << " G [" << i << ", " << j << ", " << counter << "] " << it->first << " " << dit2->first << std::endl;
                    counter++;
                }
            }
        }
    }
#endif

    //u32 tmap[buckets];
    Relation1Map *T = new Relation1Map[buckets];


    u32 **tmap_sub_bucket = new u32*[buckets];
    memset(tmap_sub_bucket, 0, sizeof(u32*) * buckets);
    for (u32 i = 0; i < buckets; i++)
    {
        tmap_sub_bucket[i] = new u32[subbuckets_T[i]];
        memset(tmap_sub_bucket[i], 0, sizeof(u32) * subbuckets_T[i]);
    }

    for (u32 i = 0; i < COL_COUNT * T_hash_entry_count; i = i + 2)
    {
        uint64_t bucket_id = hash_function(T_hashed_data[i+1]) % buckets;
        uint64_t sub_bucket_id = hash_function(T_hashed_data[i]) % subbuckets_T[bucket_id];
        tmap_bucket[bucket_id] = 1;

        auto it = T[bucket_id].find(T_hashed_data[i]);
        if( it != T[bucket_id].end() ) {
            auto it2 = (it->second)->find(T_hashed_data[i + 1]);
            if( it2 != (it->second)->end() ) {
                ;
            }
            else{
                (it->second)->insert(std::make_pair(T_hashed_data[i + 1], 0));
                T[bucket_id][T_hashed_data[i]] = it->second;
                tmap_sub_bucket[bucket_id][sub_bucket_id]++;
            }
        }
        else {
            Relation0Map *k = new Relation0Map;
            k->insert(std::make_pair(T_hashed_data[i + 1], 0));
            T[bucket_id].insert(std::make_pair(T_hashed_data[i],k));
            tmap_sub_bucket[bucket_id][sub_bucket_id]++;
        }
    }

#if 0
    for (u32 i = 0; i < buckets; i++)
    {
        for (u32 j = 0; j < subbuckets_T[i]; j++)
        {
            int counter = 0;
            for (auto it = T[i][j].begin(); it != T[i][j].end(); it++)
            {
                Relation0Map* it2 = it->second;
                for (auto dit2 = it2->begin(); dit2 != it2->end(); dit2++)
                {

                    //std::cout << "Rank " << rank << " T [" << i << ", " << j << ", " << counter << "] " << it->first << " " << dit2->first << std::endl;
                    counter++;
                }
            }
        }
    }
#endif



    Relation1Map *dT = new Relation1Map[buckets];
    memset(dT, 0, sizeof(Relation1Map) * buckets);

    int running_t_count = 0;
    for (u32 i = 0; i < COL_COUNT * T_hash_entry_count; i = i + 2)
    {
        uint64_t bucket_id = hash_function(T_hashed_data[i+1]) % buckets;

        auto it = dT[bucket_id].find(T_hashed_data[i]);
        if( it != dT[bucket_id].end() ) {
            auto it2 = (it->second)->find(T_hashed_data[i + 1]);
            if( it2 != (it->second)->end() ) {
                ;
            }
            else{
                (it->second)->insert(std::make_pair(T_hashed_data[i + 1], 0));
                dT[bucket_id].insert(std::make_pair(T_hashed_data[i], it->second));
                running_t_count++;
            }
        }
        else {
            Relation0Map *k = new Relation0Map;
            k->insert(std::make_pair(T_hashed_data[i + 1], 0));
            dT[bucket_id].insert(std::make_pair(T_hashed_data[i],k));
            running_t_count++;
        }
    }
    delete[] T_hashed_data;
    delete[] G_hashed_data;

    double relation_end = MPI_Wtime();

    double join_end = 0, join_start = 0;
    double load_balance_T_start = 0, load_balance_T_end = 0;
    double load_balance_G_start = 0, load_balance_G_end = 0;

    double lb_factor = atof(argv[6]);

    // No load balancing
    if (atoi(argv[5]) == 0)
    {
        join_start = MPI_Wtime();
        u32 g_t_bucket_count = 0;
        for (u32 bk = 0; bk < buckets; bk++)
        {
            if (gmap_bucket[bk] == 1 || tmap_bucket[bk] == 1)
                g_t_bucket_count++;
        }
        u32* g_t_bucket_indices = new u32[g_t_bucket_count];

        u32 g_t_count = 0;
        for (u32 bk = 0; bk < buckets; bk++)
        {
            if (gmap_bucket[bk] == 1 || tmap_bucket[bk] == 1)
            {
                g_t_bucket_indices[g_t_count] = bk;
                g_t_count++;
            }
        }

        int lb = 0;
        double time = 0;

        parallel_map_join(dT, G, gmap_bucket, T, tmap_bucket, 0, &lb, &running_t_count, &time, buckets, subbuckets_T, g_t_bucket_indices, g_t_bucket_count, gmap_distinct_sub_bucket_rank_count, gmap_distinct_sub_bucket_rank, tmap_distinct_sub_bucket_rank_count, tmap_distinct_sub_bucket_rank, tmap_sub_bucket_rank, tmap_sub_bucket);

        int lc = 1;
        while(true)
        {
            parallel_map_join(dT, G, gmap_bucket, T, tmap_bucket, lc, &lb, &running_t_count, &time, buckets, subbuckets_T, g_t_bucket_indices, g_t_bucket_count, gmap_distinct_sub_bucket_rank_count, gmap_distinct_sub_bucket_rank, tmap_distinct_sub_bucket_rank_count, tmap_distinct_sub_bucket_rank,tmap_sub_bucket_rank, tmap_sub_bucket);

            if (lb == 1)  break;
            lc++;
        }
        delete[] g_t_bucket_indices;
        join_end = MPI_Wtime();
    }
    // Initial load balancing
    else if (atoi(argv[5]) == 1)
    {
        if (rank == 0)
            std::cout << "[Before] Load balancing G" << std::endl;

        load_balance_G_start = MPI_Wtime();
        load_balance_G(buckets, gmap_bucket, gmap_sub_bucket, subbuckets_G, G, gmap_sub_bucket_rank, gmap_distinct_sub_bucket_rank_count, gmap_distinct_sub_bucket_rank, atoi(argv[4]), lb_factor);
        load_balance_G_end = MPI_Wtime();

        if (rank == 0)
            std::cout << "[After] Load balancing G" << std::endl;

        if (rank == 0)
            std::cout << "[Before] Load balancing T" << std::endl;

        load_balance_T_start = MPI_Wtime();
        load_balance_T(buckets, tmap_bucket, tmap_sub_bucket, subbuckets_T, T, dT, tmap_sub_bucket_rank, tmap_distinct_sub_bucket_rank_count, tmap_distinct_sub_bucket_rank, atoi(argv[4]), lb_factor);
        load_balance_T_end = MPI_Wtime();

        join_start = MPI_Wtime();

        u32 g_t_bucket_count = 0;
        for (u32 bk = 0; bk < buckets; bk++)
        {
            if (gmap_bucket[bk] == 1 || tmap_bucket[bk] == 1)
                g_t_bucket_count++;
        }
        u32* g_t_bucket_indices = new u32[g_t_bucket_count];

        u32 g_t_count = 0;
        for (u32 bk = 0; bk < buckets; bk++)
        {
            if (gmap_bucket[bk] == 1 || tmap_bucket[bk] == 1)
            {
                g_t_bucket_indices[g_t_count] = bk;
                g_t_count++;
            }
        }

        if (rank == 0)
            std::cout << "[After] Load balancing T" << std::endl;



        int lb = 0;
        double time = 0;

        parallel_map_join(dT, G, gmap_bucket, T, tmap_bucket, 0, &lb, &running_t_count, &time, buckets, subbuckets_T, g_t_bucket_indices, g_t_bucket_count, gmap_distinct_sub_bucket_rank_count, gmap_distinct_sub_bucket_rank, tmap_distinct_sub_bucket_rank_count, tmap_distinct_sub_bucket_rank, tmap_sub_bucket_rank, tmap_sub_bucket);

        int lc = 1;
        while(true)
        {
            parallel_map_join(dT, G, gmap_bucket, T, tmap_bucket, lc, &lb, &running_t_count, &time, buckets, subbuckets_T, g_t_bucket_indices, g_t_bucket_count, gmap_distinct_sub_bucket_rank_count, gmap_distinct_sub_bucket_rank, tmap_distinct_sub_bucket_rank_count, tmap_distinct_sub_bucket_rank,tmap_sub_bucket_rank, tmap_sub_bucket);

            if (lb == 1)  break;
            lc++;
        }

        delete[] g_t_bucket_indices;
        join_end = MPI_Wtime();


    }
    // load balancing every 20 iterations
    else if (atoi(argv[5]) == 2)
    {
        if (rank == 0)
            std::cout << "[Before] Load balancing G" << std::endl;

        load_balance_G_start = MPI_Wtime();
        load_balance_G(buckets, gmap_bucket, gmap_sub_bucket, subbuckets_G, G, gmap_sub_bucket_rank, gmap_distinct_sub_bucket_rank_count, gmap_distinct_sub_bucket_rank, atoi(argv[4]), lb_factor);
        load_balance_G_end = MPI_Wtime();

        if (rank == 0)
            std::cout << "[After] Load balancing G" << std::endl;

        if (rank == 0)
            std::cout << "[Before] Load balancing T" << std::endl;

        load_balance_T_start = MPI_Wtime();
        load_balance_T(buckets, tmap_bucket, tmap_sub_bucket, subbuckets_T, T, dT, tmap_sub_bucket_rank, tmap_distinct_sub_bucket_rank_count, tmap_distinct_sub_bucket_rank, atoi(argv[4]), lb_factor);
        load_balance_T_end = MPI_Wtime();

        join_start = MPI_Wtime();

        u32 g_t_bucket_count = 0;
        for (u32 bk = 0; bk < buckets; bk++)
        {
            if (gmap_bucket[bk] == 1 || tmap_bucket[bk] == 1)
                g_t_bucket_count++;
        }
        u32* g_t_bucket_indices = new u32[g_t_bucket_count];

        u32 g_t_count = 0;
        for (u32 bk = 0; bk < buckets; bk++)
        {
            if (gmap_bucket[bk] == 1 || tmap_bucket[bk] == 1)
            {
                g_t_bucket_indices[g_t_count] = bk;
                g_t_count++;
            }
        }

        if (rank == 0)
            std::cout << "[After] Load balancing T" << std::endl;

        int lb = 0;
        double time = 0;

        parallel_map_join(dT, G, gmap_bucket, T, tmap_bucket, 0, &lb, &running_t_count, &time, buckets, subbuckets_T, g_t_bucket_indices, g_t_bucket_count, gmap_distinct_sub_bucket_rank_count, gmap_distinct_sub_bucket_rank, tmap_distinct_sub_bucket_rank_count, tmap_distinct_sub_bucket_rank, tmap_sub_bucket_rank, tmap_sub_bucket);

        int lc = 1;
        while(true)
        {
            parallel_map_join(dT, G, gmap_bucket, T, tmap_bucket, lc, &lb, &running_t_count, &time, buckets, subbuckets_T, g_t_bucket_indices, g_t_bucket_count, gmap_distinct_sub_bucket_rank_count, gmap_distinct_sub_bucket_rank, tmap_distinct_sub_bucket_rank_count, tmap_distinct_sub_bucket_rank,tmap_sub_bucket_rank, tmap_sub_bucket);

            if (lc % 100 == 0)
            {
                load_balance_T(buckets, tmap_bucket, tmap_sub_bucket, subbuckets_T, T, dT, tmap_sub_bucket_rank, tmap_distinct_sub_bucket_rank_count, tmap_distinct_sub_bucket_rank, atoi(argv[4]), lb_factor);

                g_t_bucket_count = 0;
                for (u32 bk = 0; bk < buckets; bk++)
                {
                    if (gmap_bucket[bk] == 1 || tmap_bucket[bk] == 1)
                        g_t_bucket_count++;
                }

                delete[] g_t_bucket_indices;
                g_t_bucket_indices = new u32[g_t_bucket_count];

                g_t_count = 0;
                for (u32 bk = 0; bk < buckets; bk++)
                {
                    if (gmap_bucket[bk] == 1 || tmap_bucket[bk] == 1)
                    {
                        g_t_bucket_indices[g_t_count] = bk;
                        g_t_count++;
                    }
                }
            }

            if (lb == 1)  break;
            lc++;
        }

        delete[] g_t_bucket_indices;
        join_end = MPI_Wtime();
    }



    u64 total_sum = 0;
    u64 Tcounter = 0;
    for (u32 i = 0; i < buckets; i++)
    {

        for ( auto local_it = T[i].begin(); local_it!= T[i].end(); ++local_it )
        {
            Relation0Map* k = local_it->second;
            for (auto it2 = k->begin(); it2 != k->end(); it2++)
                Tcounter++;
        }

    }
    MPI_Allreduce(&Tcounter, &total_sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);

    if (rank == 0)
    {
        std::cout << "n: "
                  << nprocs
                  << " G: "
                  << global_row_count
                  << " A2a " << atoi(argv[4])
                << " T: " << total_sum
                << " Read time: " << (ior_end - ior_start)
                << " Init hash: " << (hash_end - hash_start)
                << " Relation init: " << (relation_end - relation_start)
                << " Load balG: " << (load_balance_G_end - load_balance_G_start)
                << " Load bal T: " << (load_balance_T_end - load_balance_T_start)
                << " Fixed point: " << (join_end - join_start)
                << " Total: " << (join_end - ior_start)
                << " [" << (ior_end - ior_start) + (hash_end - hash_start) + (relation_end - relation_start) + (load_balance_G_end - load_balance_G_start) + (load_balance_T_end - load_balance_T_start) + (join_end - join_start) << "]" << std::endl;
    }

    if (rank == 0)
    {
        std::cout << "n: " << nprocs
                  << " G: " << global_row_count
                  << " T: " << total_sum
                  << " CT: " << atoi(argv[4])
                << " BT: " << atoi(argv[5])
                << " Read time: " << (ior_end - ior_start)
                << " Init hash: " << (hash_end - hash_start)
                << " Relation init: " << (relation_end - relation_start)
                << " Fixed point: " << (join_end - join_start)
                << " Total: " << (join_end - ior_start)
                << " [" << (ior_end - ior_start) + (hash_end - hash_start) + (relation_end - relation_start) + (load_balance_G_end - load_balance_G_start) + (load_balance_T_end - load_balance_T_start) + (join_end - join_start) << "]" << std::endl;
    }


    for (u32 i = 0; i < buckets; i++)
    {
        Relation1Map::iterator iy1 = T[i].begin();
        for(; iy1 != T[i].end(); iy1++)
            delete (iy1->second);


        Relation1Map::iterator iy2 = G[i].begin();
        for(; iy2 != G[i].end(); iy2++)
            delete (iy2->second);


        Relation1Map::iterator iy = dT[i].begin();
        for(; iy != dT[i].end(); iy++)
            delete (iy->second);


        delete[] tmap_sub_bucket_rank[i];
        delete[] gmap_sub_bucket_rank[i];

        delete[] gmap_sub_bucket[i];
        delete[] tmap_sub_bucket[i];
    }

    delete[] gmap_distinct_sub_bucket_rank_count;
    for (u64 b = 0; b < buckets; b++)
        delete[] gmap_distinct_sub_bucket_rank[b];
    delete[] gmap_distinct_sub_bucket_rank;

    delete[] tmap_distinct_sub_bucket_rank_count;
    for (u64 b = 0; b < buckets; b++)
        delete[] tmap_distinct_sub_bucket_rank[b];
    delete[] tmap_distinct_sub_bucket_rank;


    delete[] tmap_sub_bucket_rank;
    delete[] gmap_sub_bucket_rank;

    delete[] gmap_sub_bucket;
    delete[] tmap_sub_bucket;

    delete[] dT;
    delete[] T;
    delete[] G;

    delete[] tmap_bucket;
    delete[] gmap_bucket;

    delete[] subbuckets_G;
    delete[] subbuckets_T;



#endif
    // Finalizing MPI
    MPI_Finalize();

    return 0;
}
