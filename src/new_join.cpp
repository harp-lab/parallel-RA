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
#include <vector>

#include "vector.h"

#include <unordered_set>

std::vector<u64> delta_vector;

#define DEBUG 1

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
            printf("1 Wrong input format (Meta Data)\n");
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


    *read_buffer = new u64[*local_row_count * COL_COUNT];
    u64 rb_size = pread(fp, *read_buffer, *local_row_count * COL_COUNT * sizeof(u64), read_offset * COL_COUNT * sizeof(u64));
    if (rb_size != *local_row_count * COL_COUNT * sizeof(u64))
    {
        printf("2 Wrong input format (Meta Data) %d %d\n", (int)rb_size, (int)(*local_row_count * COL_COUNT * sizeof(u64)));
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    close(fp);


    return;
}

void buffer_data_to_hash_buffer(u32 local_number_of_rows, u64* input_data,  int hash_column_index, u64** outer_hash_data, u32* outer_hash_buffer_size, MPI_Comm comm, int buckets, u32* subbuckets_G, u32* subbuckets_T, int** gmap_sub_bucket_rank)
{
    /* process_size[j] stores the number of samples to be sent to process with rank j */
    int process_size[nprocs];
    memset(process_size, 0, nprocs * sizeof(int));

    /* vector[i] contains the data that needs to be sent to process i */
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

            int index = gmap_sub_bucket_rank[bucket_id][sub_bucket_id];
            process_size[index] = process_size[index] + COL_COUNT;

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

            int index = gmap_sub_bucket_rank[bucket_id][sub_bucket_id];
            process_size[index] = process_size[index] + COL_COUNT;

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

#if 0
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



void parallel_map_join(Relation1Map*& delT, Relation1Map*& G, u32* gmap_bucket, Relation1Map*& T, u32* tmap, int lc, int* lb, u32 buckets, u32* subbuckets_T, u32* g_t_bucket_indices, u32 g_t_actual_bucket_count, int* gmap_distinct_sub_bucket_rank_count,  int** gmap_distinct_sub_bucket_rank, int* tmap_distinct_sub_bucket_rank_count, int** tmap_distinct_sub_bucket_rank, int** tmap_sub_bucket_rank, u32**& tmap_sub_bucket)
{
    double c = 0;
    double jo = 0;

    double t1_s = MPI_Wtime();

    u64 val[2] = {0, 0};

    vector_buffer* process_data_vector = (vector_buffer*)malloc(sizeof(vector_buffer) * nprocs);
    for (u32 i = 0; i < nprocs; ++i) {
        process_data_vector[i] = vector_buffer_create_empty();
    }


#if 1
    u32 local_join_duplicates = 0;
    u32 local_join_inserts = 0;
    int process_size[nprocs];
    memset(process_size, 0, nprocs * sizeof(int));
    //for (u32 i3 = 0; i3 < g_t_actual_bucket_count; i3++)
    for (u32 i = 0; i < buckets; i++)
    {
        double c_s = MPI_Wtime();
        //u32 i = g_t_bucket_indices[i3];

#if DEBUG
        if (rank == 0)
            std::cout << "[" << i << "] Meta Comm Start" << std::endl;
#endif

        vector_buffer delT_temp = vector_buffer_create_empty();

        for (auto it = T[i].begin(); it != T[i].end(); it++)
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
        std::cout << "total_buffer_size " << total_buffer_size << std::endl;
        for (u32 k1 = 0; k1 < total_buffer_size; k1=k1+2)
        {
            auto itd = G[i].find(recvbuf[k1]);
            if( itd != G[i].end() ) {
                Relation0Map* Git = itd->second;
                for (auto it2 = Git->begin(); it2 != Git->end(); it2++)
                {
                    auto itx = tempT.find(recvbuf[k1+1]);
                    if( itx != tempT.end() )
                    {
                        auto it2x = (itx->second)->find(it2->first);
                        if( it2x != (itx->second)->end() )
                        {
                            local_join_duplicates++;
                        }
                        else
                        {
                            (itx->second)->insert(std::make_pair(it2->first, 0));
                            tempT[recvbuf[k1+1]] = itx->second;
                            local_join_inserts++;

                            /* because we are sending to T and dt and they are hashed on the second column */
                            uint64_t bucket_id = hash_function(it2->first) % buckets;
                            uint64_t sub_bucket_id = hash_function(recvbuf[k1+1]) % subbuckets_T[bucket_id];

                            int index = tmap_sub_bucket_rank[bucket_id][sub_bucket_id];
                            /* uint64_t index = (hash_function((bucket_id << 32)^sub_bucket_id))%nprocs; */

                            val[0] = recvbuf[k1+1];
                            val[1] = it2->first;
                            vector_buffer_append(&process_data_vector[index], (unsigned char *) val, sizeof(u64)*2);

                            process_size[index] = process_size[index] + COL_COUNT;
                        }
                    }
                    else
                    {
                        Relation0Map* k = new Relation0Map();
                        k->insert(std::make_pair(it2->first, 0));
                        tempT[recvbuf[k1+1]] = k;
                        local_join_inserts++;

                        /* because we are sending to T and dt and they are hashed on the second column */
                        uint64_t bucket_id = hash_function(it2->first) % buckets;
                        uint64_t sub_bucket_id = hash_function(recvbuf[k1+1]) % subbuckets_T[bucket_id];

                        /*uint64_t index = (hash_function((bucket_id << 32)^sub_bucket_id))%nprocs; */
                        int index = tmap_sub_bucket_rank[bucket_id][sub_bucket_id];

                        val[0] = recvbuf[k1+1];
                        val[1] = it2->first;
                        vector_buffer_append(&process_data_vector[index], (unsigned char *) val, sizeof(u64)*2);

                        process_size[index] = process_size[index] + COL_COUNT;
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



#if 1
#if 1
    double c1 = MPI_Wtime();

#if DEBUG
    if (rank == 0)
        std::cout << "Global Comm Start" << std::endl;
#endif

    int prefix_sum_process_size[nprocs];
    memset(prefix_sum_process_size, 0, nprocs * sizeof(int));



    int process_data_buffer_size = prefix_sum_process_size[nprocs - 1] + process_size[nprocs - 1];

    u64* process_data = 0;
    process_data = new u64[process_data_buffer_size];
    memset(process_data, 0, process_data_buffer_size * sizeof(u64));

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

    /*
    for (u32 i3 = 0; i3 < g_t_actual_bucket_count; i3++)
    {
        u32 i = g_t_bucket_indices[i3];

        Relation1Map::iterator iy = delT[i].begin();
        for(; iy != delT[i].end(); iy++)
            delete (iy->second);
        delT[i].clear();
    }
    */

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

    u64 dtduplicates = 0;
    unsigned long long dtcount = 0;

    std::cout << "XXXXXXX " << outer_hash_buffer_size << std::endl;
    for (u32 i = 0; i < outer_hash_buffer_size; i = i + 2)
    {
        uint64_t bucket_id = hash_function(hash_buffer[i]) % buckets;
        uint64_t sub_bucket_id = hash_function(hash_buffer[i+1]) % subbuckets_T[bucket_id];

        auto it = delT[bucket_id].find(hash_buffer[i]);
        if ( it != delT[bucket_id].end() ) {
            auto it2 = (it->second)->find(hash_buffer[i + 1]);
            if ( it2 != (it->second)->end() ) {
                tduplicates++;
            }
            else {
                (it->second)->insert(std::make_pair(hash_buffer[i + 1], 0));
                delT[bucket_id].insert(std::make_pair(hash_buffer[i], it->second));
                tmap[bucket_id] = 1;
                tmap_sub_bucket[bucket_id][sub_bucket_id]++;

                tcount++;

            }
        }
        else {
            Relation0Map *k = new Relation0Map;
            k->insert(std::make_pair(hash_buffer[i + 1], 0));
            delT[bucket_id].insert(std::make_pair(hash_buffer[i],k));
            tmap[bucket_id] = 1;
            tmap_sub_bucket[bucket_id][sub_bucket_id]++;
            tcount++;
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
        unsigned long long sum = 0;
        MPI_Allreduce(&dtcount, &sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_BOR, comm);
        if(sum == 0)
            *lb = 1;
        else
            *lb = 0;
    }


    double v2 = MPI_Wtime();


    if (rank == 0)
    {
        std::cout << "Iteration Number: " << lc << std::endl;
        std::cout << "  Clique comm time: " << c << std::endl;

        std::cout << "  Local Join: " << jo << std::endl;
        std::cout << "    LJ (inserts) " << local_join_inserts
                  << "    LJ (duplicates) " << local_join_duplicates << std::endl;

        std::cout << "  All to All Comm: " << global_comm_time << std::endl;

        std::cout << "  Insert time: " << insert_time << std::endl;
        std::cout << "    T (inserts) " << tcount
                  << "    T (duplicates) " << tduplicates;
        std::cout << "    dT (inserts) " << dtcount
                  << "    dT (duplicates) " << dtduplicates << std::endl;

        std::cout << "  Verification time: " << (v2 - v1) << std::endl;

        std::cout << "Total time: " << c + jo + global_comm_time + insert_time + (v2-v1) << " " << (v2 - t1_s) << std::endl;
    }
#endif
#endif
}


void load_balance_G(u32 buckets, u32*& gmap_bucket, u32**& gmap_sub_bucket, u32*& subbuckets_G, Relation1Map*& G, int**& gmap_sub_bucket_rank, int*& gmap_distinct_sub_bucket_rank_count, int**& gmap_distinct_sub_bucket_rank, float lb_factor)
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
    {
        if (global_g_new_sub_bucket[i] == 0)
            global_g_new_sub_bucket[i] = 1;
#if 0
        if (global_g_new_sub_bucket[i] > (int)ceil((float)nprocs/32))
            global_g_new_sub_bucket[i] = (int)ceil((float)nprocs/32);
#endif
    }

#if DEBUG
    if (rank == 0)
        std::cout << "[G Load Balance] New sub-bucket counts" << std::endl;
#endif

    int rcount = gmap_sub_bucket_rank[buckets-1][subbuckets_G[buckets-1] - 1] + 1;

    for (u64 b = 0; b < buckets; b++)
    {
        if (subbuckets_G[b] < global_g_new_sub_bucket[b])
        {
#if 0
            gmap_bucket[b] = 0;
#endif
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
    for (int in = 0; in < outer_hash_buffer_size; in = in + 2)
    {
        uint64_t bucket_id = hash_function((outer_hash_data)[in]) % buckets;
        uint64_t sub_bucket_id = hash_function((outer_hash_data)[in+1]) % global_g_new_sub_bucket[bucket_id];
        gmap_bucket[bucket_id] = 1;

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
            if (subbuckets_G[i] != global_g_new_sub_bucket[i])
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


int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc != 6)
    {
        if (rank == 0)
        {
            std::cout << "Usage:" << "TC <filename> <bucket count factor> <initial sub bucket count> <load balancing factor>"   << std::endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }

#if 1
    comm = MPI_COMM_WORLD;
    double running_time = 0;
    double ior_start = MPI_Wtime();

    nprocs = (u32)size;

    u32 buckets;
    u32 *subbuckets_G;
    u32 *subbuckets_T;

    buckets = nprocs * atoi(argv[3]);

    subbuckets_G = new u32[buckets];
    memset(subbuckets_G, 0, buckets*sizeof(u32));

    subbuckets_T = new u32[buckets];
    memset(subbuckets_T, 0, buckets*sizeof(u32));

    for (u64 b = 0; b < buckets; b++)
    {
        subbuckets_G[b] = atoi(argv[4]);
        subbuckets_T[b] = atoi(argv[4]);
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



    u32 entry_count_G;
    u64 *input_buffer_G = NULL;
    parallel_read_input_relation_from_file_to_local_buffer(argv[1], &input_buffer_G, &entry_count_G);

    u32 entry_count_T;
    u64 *input_buffer_T = NULL;
    parallel_read_input_relation_from_file_to_local_buffer(argv[2], &input_buffer_T, &entry_count_T);

    double ior_end = MPI_Wtime();
    running_time = running_time + (ior_end - ior_start);

    if (rank == 0)
        std::cout << "Parallel IO of relations " << (ior_end - ior_start) << "[ " << running_time << " ]" <<  std::endl;

    double hash_start = MPI_Wtime();

    u32 G_hash_entry_count;
    u64 *G_hashed_data = NULL;
    buffer_data_to_hash_buffer(entry_count_G, input_buffer_G, 0, &G_hashed_data, &G_hash_entry_count, MPI_COMM_WORLD, buckets, subbuckets_G, subbuckets_T, gmap_sub_bucket_rank);

    std::cout << "[G] Rank: " << rank << " Hashed count: " << G_hash_entry_count << " Ratio: " << (float)((float) G_hash_entry_count/global_row_count)*100 << std::endl;


    u32 T_hash_entry_count;
    u64 *T_hashed_data = NULL;
    buffer_data_to_hash_buffer(entry_count_T, input_buffer_T,  1, &T_hashed_data, &T_hash_entry_count, MPI_COMM_WORLD, buckets, subbuckets_T, subbuckets_G, tmap_sub_bucket_rank);

    std::cout << "[T] Rank: " << rank << " Hashed count: " << T_hash_entry_count << " Ratio: " << (float)((float) T_hash_entry_count/global_row_count)*100 << std::endl;

    delete[] input_buffer_G;
    delete[] input_buffer_T;
    double hash_end = MPI_Wtime();

    running_time = running_time + (hash_end - hash_start);

    if (rank == 0)
        std::cout << "Initial hashing " << (hash_end - hash_start) << "[ " << running_time << " ]" <<  std::endl;

    double relation_start = MPI_Wtime();


    u32 **gmap_sub_bucket = new u32*[buckets];
    memset(gmap_sub_bucket, 0, sizeof(u32*) * buckets);

    for (u32 i = 0; i < buckets; i++)
    {
        gmap_sub_bucket[i] = new u32[subbuckets_G[i]];
        memset(gmap_sub_bucket[i], 0, sizeof(u32) * subbuckets_G[i]);
    }

    Relation1Map *G = new Relation1Map[buckets];


    u32 gcount = 0;
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
                gcount++;
            }
        }
        else {
            Relation0Map *k = new Relation0Map;
            k->insert(std::make_pair(G_hashed_data[i + 1], 0));
            G[bucket_id].insert(std::make_pair(G_hashed_data[i],k));
            gmap_sub_bucket[bucket_id][sub_bucket_id]++;
            gcount++;
        }
    }


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
        uint64_t bucket_id = hash_function(T_hashed_data[i]) % buckets;
        uint64_t sub_bucket_id = hash_function(T_hashed_data[i+1]) % subbuckets_T[bucket_id];
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

    delete[] T_hashed_data;
    delete[] G_hashed_data;


    double load_balance_T_start = 0, load_balance_T_end = 0;
    double load_balance_G_start = 0, load_balance_G_end = 0;
    //double lb_factor = atof(argv[5]);

    if (rank == 0)
        std::cout << "[Before] Load balancing G" << std::endl;

    load_balance_G_start = MPI_Wtime();
    //load_balance_G(buckets, gmap_bucket, gmap_sub_bucket, subbuckets_G, G, gmap_sub_bucket_rank, gmap_distinct_sub_bucket_rank_count, gmap_distinct_sub_bucket_rank, lb_factor);
    load_balance_G_end = MPI_Wtime();

    if (rank == 0)
        std::cout << "[After] Load balancing G" << std::endl;

    if (rank == 0)
        std::cout << "[Before] Load balancing T" << std::endl;

    load_balance_T_start = MPI_Wtime();
    //load_balance_G(buckets, tmap_bucket, tmap_sub_bucket, subbuckets_T, T, tmap_sub_bucket_rank, tmap_distinct_sub_bucket_rank_count, tmap_distinct_sub_bucket_rank, lb_factor);
    load_balance_T_end = MPI_Wtime();



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
            std::cout << "XXXXXXXX rank " << rank << " " << bk << std::endl;
            g_t_count++;
        }
    }

    double relation_end = MPI_Wtime();


    double join_end = 0, join_start = 0;
    join_start = MPI_Wtime();
    int lb = 0;


    int lc = 0;

    std::cout << "buckets " << buckets << std::endl;

    Relation1Map *dT = new Relation1Map[buckets];
    parallel_map_join(dT, G, gmap_bucket, T, tmap_bucket, lc, &lb, buckets, subbuckets_T, g_t_bucket_indices, g_t_bucket_count, gmap_distinct_sub_bucket_rank_count, gmap_distinct_sub_bucket_rank, tmap_distinct_sub_bucket_rank_count, tmap_distinct_sub_bucket_rank,tmap_sub_bucket_rank, tmap_sub_bucket);




    delete[] g_t_bucket_indices;
    join_end = MPI_Wtime();


    u64 total_sum = 0;
    u64 Tcounter = 0;
    for (u32 i = 0; i < buckets; i++)
    {

        for ( auto local_it = dT[i].begin(); local_it!= dT[i].end(); ++local_it )
        {
            Relation0Map* k = local_it->second;
            for (auto it2 = k->begin(); it2 != k->end(); it2++)
                Tcounter++;
        }

    }
    MPI_Allreduce(&Tcounter, &total_sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);

#if 1
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
                  << " Load balG: " << (load_balance_G_end - load_balance_G_start)
                  << " Load bal T: " << (load_balance_T_end - load_balance_T_start)
                  << " Fixed point: " << (join_end - join_start)
                  << " Total: " << (join_end - ior_start)
                  << " [" << (ior_end - ior_start) + (hash_end - hash_start) + (relation_end - relation_start) + (load_balance_G_end - load_balance_G_start) + (load_balance_T_end - load_balance_T_start) + (join_end - join_start) << "]" << std::endl;
    }
#endif

    if (rank == 0)
    {
        std::cout << "n: " << nprocs
                  << " G: " << global_row_count
                  << " T: " << total_sum
                  << " Read time: " << (ior_end - ior_start)
                  << " Init hash: " << (hash_end - hash_start)
                  << " Relation init: " << (relation_end - relation_start)
                  << " Fixed point: " << (join_end - join_start)
                  << " Total: " << (join_end - ior_start)
                  << " [" << (ior_end - ior_start) + (hash_end - hash_start) + (relation_end - relation_start) + (join_end - join_start) << "]" << std::endl;
    }


    if (rank == 0)
    {
#if 0
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
#endif
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
    MPI_Finalize();

    return 0;
}
