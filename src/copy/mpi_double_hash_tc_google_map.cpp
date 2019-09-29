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

//#include "btree.h"
//#include "btree_relation.h"
#include "btree/btree_map.h"

typedef btree::btree_map<u64, u64> Relation0Map;
typedef btree::btree_map<u64, btree::btree_map<u64, u64>* > Relation1Map;

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
    u64 rb_size = pread(fp, *read_buffer, *local_row_count * global_col_count * sizeof(u64), read_offset * global_col_count * sizeof(u64));
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
    u64 buckets = ceil((float)(3 * nprocs) / 4);
    u64 subbuckets = 2;

    std::cout << "Number of processes: " << nprocs << " Bucket count: " << buckets << std::endl;

    /* process_size[j] stores the number of samples to be sent to process with rank j */
    int process_size[nprocs];
    memset(process_size, 0, nprocs * sizeof(int));

    /* vector[i] contains the data that needs to be sent to process i */
    std::vector<u64> process_data_vector[nprocs];

    if (hash_column_index == 0)
    {
        for (u32 i = 0; i < local_number_of_rows * col_count; i=i+2)
        {
            uint64_t bucket_id = outer_hash(input_data[i]) % buckets;
            uint64_t sub_bucket_id = outer_hash(input_data[i+1]) % subbuckets;

            uint64_t index = outer_hash(bucket_id*sub_bucket_id)%nprocs;
            process_size[index] = process_size[index] + col_count;

            process_data_vector[index].push_back(input_data[i]);
            process_data_vector[index].push_back(input_data[i + 1]);
        }
    }
    else
    {
        for (u32 i = 0; i < local_number_of_rows * col_count; i=i+2)
        {
            uint64_t bucket_id = outer_hash(input_data[i]) % buckets;
            uint64_t sub_bucket_id = outer_hash(input_data[i+1]) % subbuckets;

            uint64_t index = outer_hash(bucket_id*sub_bucket_id)%nprocs;
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

    return;
}


Relation1Map *** parallel_map_join(Relation1Map*** delT, u32* dtmap, Relation1Map**& G, u32* gmap, Relation1Map**& T, int lc, int* lb, int* running_t_count, double* running_time)
{

    u64 buckets = ceil((float)(3 * nprocs) / 4);
    u64 subbuckets = 2;


    Relation1Map ***delTT = new Relation1Map**[buckets];
    for (u32 i = 0; i < buckets; i++)
    {
        delTT[i] = new Relation1Map*[subbuckets];
        for (u32 j = 0; j < subbuckets; j++)
            delTT[i][j] = new Relation1Map;
    }

    double j1 = 0, j2 = 0;
    double c1 = 0, c2 = 0;
    double i1 = 0, i2 = 0;
    double v1 = 0, v2 = 0;

    //j1 = MPI_Wtime();

    tuple<2> t;
    t[0] = -1;
    t[1] = -1;
    //tuple<2> selectall(t);
    Relation1Map tempT;

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

    tuple<2> dt;
    tuple<2> s;
    s[1] = -1;
    //tuple<2> select;
    //uint64_t index;

    //tuple<2> dt;
    //u64 k;

    std::vector<u64> delT_temp[buckets];
    int color[buckets];
    memset(color, 0, buckets * sizeof(int));

    MPI_Comm group_comm[nprocs];

    for (u32 i = 0; i < buckets; i++)
    {
        if (dtmap[i] == 1 || gmap[i] == 1)
            color[i] = 1;

        for (u32 j = 0; j < subbuckets; j++)
        {
            for (auto it = delT[i][j]->begin(); it != delT[i][j]->end(); it++)
            {
                Relation0Map* it2 = it->second;
                for (auto dit2 = it2->begin(); dit2 != it2->end(); dit2++)
                {
                    delT_temp[i].push_back(it->first);
                    delT_temp[i].push_back(dit2->first);

                    //std::cout << "TEST: " << it->first << " " << dit2->first << std::endl;
                }
            }
        }

        MPI_Comm_split(MPI_COMM_WORLD, color[i], rank, &(group_comm[i]));
    }


    int count = 0;
    int gnprocs = 0;
    u64 tduplicates = 0;
    for (u32 i = 0; i < buckets; i++)
    {
        MPI_Comm_size(group_comm[i], &gnprocs);

        int process_size[gnprocs];
        memset(process_size, 0, gnprocs * sizeof(int));

        int recvcounts[gnprocs];
        int displs[gnprocs];
        for (int n = 0; n < gnprocs; n++)
        {
            recvcounts[n] = 1;
            displs[n] = n;
        }

        /* This step prepares for actual data transfer */
        /* Every process sends to every other process the amount of data it is going to send */
        int recv_process_size_buffer_local[gnprocs];
        memset(recv_process_size_buffer_local, 0, gnprocs * sizeof(int));

        u32 buffer_size = delT_temp[i].size();
        MPI_Allgatherv(&buffer_size, 1, MPI_INT, recv_process_size_buffer_local, recvcounts, displs, MPI_INT, group_comm[i]);

        int recv_process_prefix[gnprocs];
        memset(recv_process_prefix, 0, gnprocs * sizeof(int));

        u64 total_buffer_size = 0;
        for (int n = 0; n < gnprocs; n++)
        {
            recv_process_prefix[n] = total_buffer_size;
            total_buffer_size = total_buffer_size + recv_process_size_buffer_local[n];
        }

        u64 *recvbuf = new u64[total_buffer_size];
        MPI_Allgatherv(&delT_temp[i][0], delT_temp[i].size(), MPI_UNSIGNED_LONG_LONG, recvbuf, recv_process_size_buffer_local, recv_process_prefix, MPI_UNSIGNED_LONG_LONG, group_comm[i]);

        for (u32 j = 0; j < subbuckets; j++)
        {
            for (u32 k1 = 0; k1 < total_buffer_size; k1=k1+2)
            {
                auto itd = G[i][j].find(recvbuf[k1+1]);
                if( itd != G[i][j].end() ) {
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

                                dt[0] = recvbuf[k1];
                                dt[1] = it2->first;

                                uint64_t bucket_id = outer_hash(dt[0]) % buckets;
                                uint64_t sub_bucket_id = outer_hash(dt[1]) % subbuckets;
                                uint64_t index = outer_hash(bucket_id*sub_bucket_id)%nprocs;
                                process_data_vector[index].push_back(dt);

                                //std::cout << "I1 " << dt[0] << " " << dt[1] << std::endl;
                            }
                        }
                        else {
                            Relation0Map* k = new Relation0Map();
                            k->insert(std::make_pair(it2->first, 0));
                            //tempT.insert(std::make_pair(it->first,k));
                            tempT[recvbuf[k1]] = k;

                            dt[0] = recvbuf[k1];
                            dt[1] = it2->first;

                            uint64_t bucket_id = outer_hash(dt[0]) % buckets;
                            uint64_t sub_bucket_id = outer_hash(dt[1]) % subbuckets;
                            uint64_t index = outer_hash(bucket_id*sub_bucket_id)%nprocs;
                            process_data_vector[index].push_back(dt);

                            //std::cout << "I2 " << dt[0] << " " << dt[1] << std::endl;
                        }
                    }
                }
            }
        }

        Relation1Map::iterator ix = tempT.begin();
        for(; ix != tempT.end(); ix++)
            delete (ix->second);
        //j2 = MPI_Wtime();

        //c1 = MPI_Wtime();

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


        //c2 = MPI_Wtime();

        //i1 = MPI_Wtime();
#if 1
        for (u32 i = 0; i < outer_hash_buffer_size; i = i + 2)
        {
            uint64_t bucket_id = outer_hash(hash_buffer[i]) % buckets;
            uint64_t sub_bucket_id = outer_hash(hash_buffer[i+1]) % subbuckets;

            auto it = T[bucket_id][sub_bucket_id].find(hash_buffer[i]);
            if ( it != T[bucket_id][sub_bucket_id].end() ) {
                auto it2 = (it->second)->find(hash_buffer[i + 1]);
                if ( it2 != (it->second)->end() ) {
                    tduplicates++;
                }
                else{
                    (it->second)->insert(std::make_pair(hash_buffer[i + 1], 0));
                    T[bucket_id][sub_bucket_id].insert(std::make_pair(hash_buffer[i], it->second));

                    std::cout << "IJ 1 " << bucket_id << " " << sub_bucket_id << " " << hash_buffer[i] << " " << hash_buffer[i + 1] << std::endl;
                    tcount++;
                    {
                    auto itx = delTT[bucket_id][sub_bucket_id]->find(hash_buffer[i]);
                    if ( itx != delTT[bucket_id][sub_bucket_id]->end() ) {
                        auto it2x = (itx->second)->find(hash_buffer[i + 1]);
                        if ( it2x != (itx->second)->end() ) {
                            ;
                        }
                        else{
                            (itx->second)->insert(std::make_pair(hash_buffer[i + 1], 0));
                            delTT[bucket_id][sub_bucket_id]->insert(std::make_pair(hash_buffer[i], itx->second));
                            count++;

                            std::cout << "IJD 1 " << bucket_id << " " << sub_bucket_id << " " << hash_buffer[i] << " " << hash_buffer[i + 1] << std::endl;
                        }
                    }
                    else {
                        Relation0Map *k = new Relation0Map;
                        k->insert(std::make_pair(hash_buffer[i + 1], 0));
                        delTT[bucket_id][sub_bucket_id]->insert(std::make_pair(hash_buffer[i],k));
                        count++;

                        std::cout << "IJD 2 " << bucket_id << " " << sub_bucket_id << " " << hash_buffer[i] << " " << hash_buffer[i + 1] << std::endl;
                    }
                    }
                }
            }
            else {
                Relation0Map *k = new Relation0Map;
                k->insert(std::make_pair(hash_buffer[i + 1], 0));
                T[bucket_id][sub_bucket_id].insert(std::make_pair(hash_buffer[i],k));
                tcount++;

                std::cout << "IJ 2 " << bucket_id << " " << sub_bucket_id << " " << hash_buffer[i] << " " << hash_buffer[i + 1] << std::endl;
                {
                auto itx = delTT[bucket_id][sub_bucket_id]->find(hash_buffer[i]);
                if ( itx != delTT[bucket_id][sub_bucket_id]->end() ) {
                    auto it2x = (itx->second)->find(hash_buffer[i + 1]);
                    if ( it2x != (itx->second)->end() ) {
                        ;
                    }
                    else{
                        (itx->second)->insert(std::make_pair(hash_buffer[i + 1], 0));
                        delTT[bucket_id][sub_bucket_id]->insert(std::make_pair(hash_buffer[i], itx->second));
                        count++;

                        std::cout << "IJD 1 " << bucket_id << " " << sub_bucket_id << " " << hash_buffer[i] << " " << hash_buffer[i + 1] << std::endl;
                    }
                }
                else {
                    Relation0Map *k = new Relation0Map;
                    k->insert(std::make_pair(hash_buffer[i + 1], 0));
                    delTT[bucket_id][sub_bucket_id]->insert(std::make_pair(hash_buffer[i],k));
                    count++;

                    std::cout << "IJD 2 " << bucket_id << " " << sub_bucket_id << " " << hash_buffer[i] << " " << hash_buffer[i + 1] << std::endl;
                }
                }
            }
        }
#endif

#if 0
        for(u32 ko = 0; ko < outer_hash_buffer_size; ko = ko + 2)
        {
            uint64_t bucket_id = outer_hash(hash_buffer[ko]) % buckets;
            uint64_t sub_bucket_id = outer_hash(hash_buffer[ko+1]) % subbuckets;

            auto it = T[bucket_id][sub_bucket_id].find(hash_buffer[ko]);

            if( it != T[bucket_id][sub_bucket_id].end() ) {
                auto it2 = (it->second)->find(hash_buffer[ko + 1]);
                if( it2 != (it->second)->end() ) {
                    tduplicates++;
                }
                else
                {
                    (it->second)->insert(std::make_pair(hash_buffer[ko + 1], 0));
                    T[bucket_id][sub_bucket_id][hash_buffer[ko]] = it->second;
                    tcount++;

                    //std::cout << "Verify: " << bucket_id << " " << sub_bucket_id << std::endl;

                    //
                    auto sit = delTT[bucket_id][sub_bucket_id]->find(hash_buffer[ko]);
                    if( sit != delTT[bucket_id][sub_bucket_id]->end() ) {
                        (sit->second)->insert(std::make_pair(hash_buffer[ko + 1], 0));
                        (*delTT)[bucket_id][sub_bucket_id][hash_buffer[ko]] = sit->second;

                        //std::cout << "1    ToI " << hash_buffer[ko] << " " << hash_buffer[ko + 1] << std::endl;
                        //std::cout << "1    ij " << bucket_id << " " << sub_bucket_id << std::endl;

                        count++;
                    }
                    else {
                        //std::cout << "2    ToI " << hash_buffer[ko] << " " << hash_buffer[ko + 1] << std::endl;
                        //std::cout << "2    ij " << bucket_id << " " << sub_bucket_id << std::endl;

                        Relation0Map *kn = new Relation0Map;
                        kn->insert(std::make_pair(hash_buffer[ko + 1], 0));
                        //(*delTT)[bucket_id][sub_bucket_id][hash_buffer[ko]] = kn;
                        delTT[bucket_id][sub_bucket_id]->insert(std::make_pair(hash_buffer[ko], kn));
                        count++;
                    }
                    //
                }
            }
            else {
                Relation0Map *k2 = new Relation0Map;
                k2->insert(std::make_pair(hash_buffer[ko + 1], 0));
                T[bucket_id][sub_bucket_id][hash_buffer[ko]] = k2;
                //tcount++;

                //
                auto sit = delTT[bucket_id][sub_bucket_id]->find(hash_buffer[ko]);
                if( sit != delTT[bucket_id][sub_bucket_id]->end() ) {
                    (sit->second)->insert(std::make_pair(hash_buffer[ko + 1], 0));
                    //(*delTT).insert(std::make_pair(hash_buffer[ko],sit->second));
                    (*delTT)[bucket_id][sub_bucket_id][hash_buffer[ko]] = sit->second;
                    count++;
                }
                else {
                    Relation0Map *k = new Relation0Map;
                    k->insert(std::make_pair(hash_buffer[ko + 1], 0));
                    //(*delTT).insert(std::make_pair(hash_buffer[ko],k));
                    (*delTT)[bucket_id][sub_bucket_id][hash_buffer[ko]] = k;
                    count++;
                }
                //
            }
        }
#endif

        //myfile.close();
        //dfile.close();

        delete[] hash_buffer;
        delete[] process_data;
        //i2 = MPI_Wtime();

    }
#if 1
    //v1 = MPI_Wtime();
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
    //

    for (u32 i = 0; i < buckets; i++)
    {
        for (u32 j = 0; j < subbuckets; j++)
        {
            Relation1Map::iterator iy = delT[i][j]->begin();
            for(; iy != delT[i][j]->end(); iy++)
                delete (iy->second);
            delete delT[i][j];
        }
    }
#endif


    return delTT;

}

//create big 2d array bucket/sub-bucket most will be empty
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

    u64 buckets = ceil((float)(3 * nprocs) / 4);
    u64 subbuckets = 2;

    u32 entry_count;
    u64 *input_buffer = NULL;
    parallel_read_input_relation_from_file_to_local_buffer(argv[1], &input_buffer, &entry_count, &col_count);

    double ior_end = MPI_Wtime();

    double hash_start = MPI_Wtime();

    u32 G_hash_entry_count;
    u64 *G_hashed_data = NULL;
    buffer_data_to_hash_buffer(entry_count, col_count, input_buffer, 0, &G_hashed_data, &G_hash_entry_count, MPI_COMM_WORLD);

    std::cout << "[G] Rank: " << rank << " Hashed count: " << G_hash_entry_count << " Ratio: " << (float)((float) G_hash_entry_count/global_row_count)*100 << std::endl;


    u32 T_hash_entry_count;
    u64 *T_hashed_data = NULL;
    buffer_data_to_hash_buffer(entry_count, col_count, input_buffer,  1, &T_hashed_data, &T_hash_entry_count, MPI_COMM_WORLD);

    std::cout << "[T] Rank: " << rank << " Hashed count: " << T_hash_entry_count << " Ratio: " << (float)((float) T_hash_entry_count/global_row_count)*100 << std::endl;

    delete[] input_buffer;
    double hash_end = MPI_Wtime();


    double relation_start = MPI_Wtime();

    u32 gmap[buckets];
    memset(gmap, 0, sizeof(u32) * buckets);

    Relation1Map **G = new Relation1Map*[buckets];
    for (u32 i = 0; i < buckets; i++)
        G[i] = new Relation1Map[subbuckets];


    for (u32 i = 0; i < col_count * G_hash_entry_count; i = i + 2)
    {
        uint64_t bucket_id = outer_hash(G_hashed_data[i]) % buckets;
        uint64_t sub_bucket_id = outer_hash(G_hashed_data[i+1]) % subbuckets;
        gmap[bucket_id] = 1;

        auto it = G[bucket_id][sub_bucket_id].find(G_hashed_data[i]);
        if( it != G[bucket_id][sub_bucket_id].end() ) {
            auto it2 = (it->second)->find(G_hashed_data[i + 1]);
            if( it2 != (it->second)->end() ) {
                ;
            }
            else{
                (it->second)->insert(std::make_pair(G_hashed_data[i + 1], 0));
                G[bucket_id][sub_bucket_id][G_hashed_data[i]] = it->second;
            }
        }
        else {
            Relation0Map *k = new Relation0Map;
            k->insert(std::make_pair(G_hashed_data[i + 1], 0));
            G[bucket_id][sub_bucket_id].insert(std::make_pair(G_hashed_data[i],k));
        }
    }


    for (u32 i = 0; i < buckets; i++)
    {
        for (u32 j = 0; j < subbuckets; j++)
        {
            int counter = 0;
            for (auto it = G[i][j].begin(); it != G[i][j].end(); it++)
            {
                Relation0Map* it2 = it->second;
                for (auto dit2 = it2->begin(); dit2 != it2->end(); dit2++)
                {
                    std::cout << "G [" << i << ", " << j << ", " << counter << "] " << it->first << " " << dit2->first << std::endl;
                    counter++;
                }
            }
        }
    }


    //u32 tmap[buckets];
    Relation1Map **T = new Relation1Map*[buckets];
    for (u32 i = 0; i < buckets; i++)
        T[i] = new Relation1Map[subbuckets];

    for (u32 i = 0; i < col_count * T_hash_entry_count; i = i + 2)
    {
        uint64_t bucket_id = outer_hash(T_hashed_data[i]) % buckets;
        uint64_t sub_bucket_id = outer_hash(T_hashed_data[i+1]) % subbuckets;
        //tmap[buckets] = 1;

        auto it = T[bucket_id][sub_bucket_id].find(T_hashed_data[i]);
        if( it != T[bucket_id][sub_bucket_id].end() ) {
            auto it2 = (it->second)->find(T_hashed_data[i + 1]);
            if( it2 != (it->second)->end() ) {
                ;
            }
            else{
                (it->second)->insert(std::make_pair(T_hashed_data[i + 1], 0));
                T[bucket_id][sub_bucket_id][T_hashed_data[i]] = it->second;
            }
        }
        else {
            Relation0Map *k = new Relation0Map;
            k->insert(std::make_pair(T_hashed_data[i + 1], 0));
            T[bucket_id][sub_bucket_id].insert(std::make_pair(T_hashed_data[i],k));
        }
    }

    for (u32 i = 0; i < buckets; i++)
    {
        for (u32 j = 0; j < subbuckets; j++)
        {
            int counter = 0;
            for (auto it = T[i][j].begin(); it != T[i][j].end(); it++)
            {
                Relation0Map* it2 = it->second;
                for (auto dit2 = it2->begin(); dit2 != it2->end(); dit2++)
                {
                    std::cout << "T [" << i << ", " << j << ", " << counter << "] " << it->first << " " << dit2->first << std::endl;
                    counter++;
                }
            }
        }
    }


    u32 dtmap[buckets];
    memset(dtmap, 0, sizeof(u32) * buckets);

    Relation1Map ***dT = new Relation1Map**[buckets];
    for (u32 i = 0; i < buckets; i++)
    {
        dT[i] = new Relation1Map*[subbuckets];
        for (u32 j = 0; j < subbuckets; j++)
            dT[i][j] = new Relation1Map;
    }

    int running_t_count = 0;

    for (u32 i = 0; i < col_count * T_hash_entry_count; i = i + 2)
    {
        uint64_t bucket_id = outer_hash(T_hashed_data[i]) % buckets;
        uint64_t sub_bucket_id = outer_hash(T_hashed_data[i+1]) % subbuckets;
        dtmap[bucket_id] = 1;

        auto it = dT[bucket_id][sub_bucket_id]->find(T_hashed_data[i]);
        if( it != dT[bucket_id][sub_bucket_id]->end() ) {
            auto it2 = (it->second)->find(T_hashed_data[i + 1]);
            if( it2 != (it->second)->end() ) {
                ;
            }
            else{
                (it->second)->insert(std::make_pair(T_hashed_data[i + 1], 0));
                (*dT)[bucket_id][sub_bucket_id][T_hashed_data[i]] = it->second;
                running_t_count++;
            }
        }
        else {
            Relation0Map *k = new Relation0Map;
            k->insert(std::make_pair(T_hashed_data[i + 1], 0));
            dT[bucket_id][sub_bucket_id]->insert(std::make_pair(T_hashed_data[i],k));
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
    dT = parallel_map_join(dT, dtmap, G, gmap, T, 0, &lb, &running_t_count, &time);

#if 1
    int lc = 1;
    while(true)
    {
        dT = parallel_map_join(dT, dtmap, G, gmap, T, lc, &lb, &running_t_count, &time);

        if (lb == 1)  break;
        lc++;
    }


    for (u32 i = 0; i < buckets; i++)
    {
        for (u32 j = 0; j < subbuckets; j++)
        {
            Relation1Map::iterator iy = dT[i][j]->begin();
            for(; iy != dT[i][j]->end(); iy++)
                delete (iy->second);
            delete dT[i][j];
        }
    }



    u64 total_sum = 0;

    /*
    u64 Tcounter = 0;
    for ( auto local_it = T.begin(); local_it!= T.end(); ++local_it )
    {
        Relation0Map* k = local_it->second;
        for (auto it2 = k->begin(); it2 != k->end(); it2++)
            Tcounter++;
    }
    MPI_Allreduce(&Tcounter, &total_sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
    */
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

    for (u32 i = 0; i < buckets; i++)
    {
        for (u32 j = 0; j < subbuckets; j++)
        {
            Relation1Map::iterator iy1 = T[i][j].begin();
            for(; iy1 != T[i][j].end(); iy1++)
                delete (iy1->second);
        }
    }

    for (u32 i = 0; i < buckets; i++)
    {
        for (u32 j = 0; j < subbuckets; j++)
        {
            Relation1Map::iterator iy2 = G[i][j].begin();
            for(; iy2 != G[i][j].end(); iy2++)
                delete (iy2->second);
        }
    }

#endif
    // Finalizing MPI
    MPI_Finalize();

    return 0;
}
