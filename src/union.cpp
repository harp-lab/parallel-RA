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


#include "btree/btree_map.h"
#include "compat.h"
#include "tuple.h"

/* #define OLD 1 */
double j1x, j2;
double c1 = 0, c2 = 0;
double in1 = 0, i2 = 0;

#if OLD
#include "btree.h"
#include "btree_relation.h"
#else
typedef btree::btree_map<u64, u64> Relation0Map;
typedef btree::btree_map<u64, btree::btree_map<u64, u64>* > Relation1Map;
#endif

#define COL_COUNT 2

static int rank = 0;
static u32 nprocs = 1;
static MPI_Comm comm;


#if 1

const uint32_t Prime = 0x01000193;
const uint32_t Seed  = 0x811C9DC5;

inline uint32_t fnv1a(unsigned char oneByte, uint32_t hash = Seed)
{
  return (oneByte ^ hash) * Prime;
}

uint32_t fnv1a(const void* data, size_t numBytes, uint32_t hash = Seed)
{
  const unsigned char* ptr = (const unsigned char*)data;
  while (numBytes--)
    hash = fnv1a(*ptr++, hash);
  return hash;
}

u64 outer_hash(const u64 val)
{
    return fnv1a((u8*)(&val),sizeof(u64));
}
#endif

#if 0
void parallel_read_input_relation_from_file_to_local_buffer(const char *file_name, u64** read_buffer, u64 row_count, u64* local_row_count)
{
    u64 read_offset;
    read_offset = ceil((float)row_count / nprocs) * rank;

    if (read_offset > row_count)
    {
        *local_row_count = 0;
        return;
    }

    if (read_offset + ceil((float)row_count / nprocs) > row_count)
        *local_row_count = row_count - read_offset;
    else
        *local_row_count = (u64) ceil((float)row_count / nprocs);

    if (*local_row_count == 0)
        return;

    char data_filename[1024];
    sprintf(data_filename, "%s/data.raw", file_name);
    int fp = open(data_filename, O_RDONLY);


    size_t bsize = (u64) *local_row_count *  (u64)COL_COUNT * sizeof(u64);
    *read_buffer = new u64[*local_row_count * COL_COUNT];


    u64 rb_size = pread(fp, *read_buffer, bsize, read_offset * COL_COUNT * sizeof(u64));
    if (rb_size != *local_row_count * COL_COUNT * sizeof(u64))
    {
        std::cout << "Wrong IO: rank: " << rank << " " << rb_size << " " <<  bsize << " " << *local_row_count << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    close(fp);


    /*
    std::ifstream myFile (data_filename, std::ios::in | std::ios::binary);
    myFile.seekg (read_offset * COL_COUNT * sizeof(u64));
    myFile.read ((char*)*read_buffer, bsize);
    myFile.close();
    */


    return;
}

void buffer_data_to_hash_buffer(u64 local_number_of_rows, u64* input_data, u64** outer_hash_data, u64* outer_hash_buffer_size, MPI_Comm comm, u32 index)
{
    /* process_size[j] stores the number of samples to be sent to process with rank j */
    int process_size[nprocs];
    memset(process_size, 0, nprocs * sizeof(int));

    /* vector[i] contains the data that needs to be sent to process i */
    std::vector<tuple<2>> *process_data_vector;
    process_data_vector = new std::vector<tuple<2>>[nprocs];

    if (rank == 0 && index == 1)
        std::cout << "Local number of rows: " << local_number_of_rows << std::endl;

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
    if (rank == 0 && index == 1)
        std::cout << "PS 0" << " " << process_size[0] << std::endl;

    for (u32 i = 1; i < nprocs; i++)
    {
        prefix_sum_process_size[i] = prefix_sum_process_size[i - 1] + process_size[i - 1];
        process_size[i] = 2 * process_data_vector[i].size();

        if (rank == 0 && index == 1)
            std::cout << "PS " << i << " " << process_size[i] << std::endl;
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
    if (rank == 0 && index == 1)
        std::cout << "OHBS 0" << " " << *outer_hash_buffer_size << std::endl;
    for (u32 i = 1; i < nprocs; i++)
    {
        prefix_sum_recv_process_size_buffer[i] = prefix_sum_recv_process_size_buffer[i - 1] + recv_process_size_buffer[i - 1];
        *outer_hash_buffer_size = *outer_hash_buffer_size + recv_process_size_buffer[i];

        if (rank == 0 && index == 1)
            std::cout << "OHBS " << i << " " << *outer_hash_buffer_size << std::endl;
    }

    *outer_hash_data = new u64[*outer_hash_buffer_size];
    MPI_Alltoallv(process_data, process_size, prefix_sum_process_size, MPI_UNSIGNED_LONG_LONG, *outer_hash_data, recv_process_size_buffer, prefix_sum_recv_process_size_buffer, MPI_UNSIGNED_LONG_LONG, comm);

    *outer_hash_buffer_size = *outer_hash_buffer_size / COL_COUNT;

    delete[] process_data;

    return;
}
#endif


#if OLD
u64 parallel_join(relation<2>& O, relation<2>& G, relation<2>& T)
#else
u64 parallel_join(Relation1Map& O, Relation1Map& G, Relation1Map& T)
#endif
{
    j1x = MPI_Wtime();

    tuple<2> t;
    t[0] = -1;
    t[1] = -1;

#if OLD
    tuple<2> selectall(t);
    relation<2> tempT;
    tuple<2> select;
    u64 tuple_count = 0;
#else
    Relation1Map tempT;
#endif

    /* process_size[j] stores the number of samples to be sent to process with rank j */
    int process_size[nprocs];
    memset(process_size, 0, nprocs * sizeof(int));

    /* vector[i] contains the data that needs to be sent to process i */
    std::vector<tuple<2>> *process_data_vector;
    process_data_vector = new std::vector<tuple<2>>[nprocs];


    u64 non_deduplicate_tuple_count = 0;

    tuple<2> dt;
    tuple<2> s;
    s[1] = -1;

    uint64_t index;
    u32 join_f = 0;
    u32 join_s = 0;

#if OLD
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
                join_s++;
            }
            else
                join_f++;
        }
    }
#else

    for (auto itm = T.begin(); itm != T.end(); itm++)
    {
        Relation0Map* it2m = itm->second;
        for (auto dit2m = it2m->begin(); dit2m != it2m->end(); dit2m++)
        {
            auto itd = G.find(itm->first);
            if( itd != G.end() ) {
                Relation0Map* Git = itd->second;
                for (auto it2 = Git->begin(); it2 != Git->end(); it2++)
                {
                    auto itx = tempT.find(dit2m->first);
                    if( itx != tempT.end() )
                    {
                        auto it2x = (itx->second)->find(it2->first);
                        if( it2x != (itx->second)->end() )
                        {
                            join_f++;
                        }
                        else
                        {
                            (itx->second)->insert(std::make_pair(it2->first, 0));
                            tempT[dit2m->first] = itx->second;

                            index = outer_hash(dit2m->first)%nprocs;

                            dt[0] = dit2m->first;
                            dt[1] = it2->first;
                            process_data_vector[index].push_back(dt);
                            join_s++;
                        }
                    }
                    else
                    {
                        Relation0Map* k = new Relation0Map();
                        k->insert(std::make_pair(it2->first, 0));
                        tempT[dit2m->first] = k;

                        index = outer_hash(dit2m->first)%nprocs;

                        dt[0] = dit2m->first;
                        dt[1] = it2->first;
                        process_data_vector[index].push_back(dt);
                        join_s++;
                    }
                }
            }
        }
    }
    Relation1Map::iterator ix = tempT.begin();
    for(; ix != tempT.end(); ix++)
        delete (ix->second);

#endif
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

    in1 = MPI_Wtime();
    u64 tcount = 0;

#if OLD
    for (u32 k = 0; k < outer_hash_buffer_size; k = k + 2)
    {
        tuple<2> dt;
        dt[0] = hash_buffer[k];
        dt[1] = hash_buffer[k + 1];

        if (O.insert(dt) == true)
            tcount++;
    }
#else
    for (u32 i = 0; i < outer_hash_buffer_size; i = i + 2)
    {
        auto it = O.find(hash_buffer[i]);
        if ( it != O.end() ) {
            auto it2 = (it->second)->find(hash_buffer[i + 1]);
            if ( it2 != (it->second)->end() ) {
            }
            else {
                (it->second)->insert(std::make_pair(hash_buffer[i + 1], 0));
                O.insert(std::make_pair(hash_buffer[i], it->second));
                tcount++;
            }
        }
        else {
            Relation0Map *k = new Relation0Map;
            k->insert(std::make_pair(hash_buffer[i + 1], 0));
            O.insert(std::make_pair(hash_buffer[i],k));
            tcount++;
        }
    }
#endif
    delete[] hash_buffer;
    delete[] process_data;
    i2 = MPI_Wtime();


    if (rank == 0)
        std::cout << "FFFFF Join: " << (j2 - j1x)
                  << " Comm: " << (c2 - c1)
                  << " Insert: " << (i2 - in1)
                  << " Total: " << (i2 - j1x)
                  << " [" << (j2 - j1x) + (c2 - c1) + (i2 - in1) << "]"
                  << " Delta: " << tcount
                  << std::endl;

    return tcount;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    comm = MPI_COMM_WORLD;

    nprocs = (u32)size;

    u32 relation_count = atoi(argv[1]);
    u64 *  row_count = new u64[relation_count];
    u64 *  hash_entry_count = new u64[relation_count];

    u64 *  local_entry_count = new u64[relation_count];
    u64 ** input_buffer = new u64*[relation_count];
    u64 ** hashed_data = new u64*[relation_count];

    double * ior_start = new double[relation_count];
    double * ior_end = new double[relation_count];
    double * hash_start = new double[relation_count];
    double * hash_end = new double[relation_count];

    double * comm_start = new double[relation_count];
    double * comm_end = new double[relation_count];

    double * union_start = new double[relation_count];
    double * union_end = new double[relation_count];


    Relation1Map unionG;
    tuple<2> t1;


    u64 total_tuple = 0;
    u64 inserted_tuple = 0;
    for (u32 i1 = 0; i1 < relation_count; i1++)
    {
        ior_start[i1] = MPI_Wtime();
        row_count[i1] = atoi(argv[3 + 2*i1]);

        u64 read_offset;
        read_offset = ceil((float)row_count[i1] / nprocs) * rank;

        if (read_offset > row_count[i1])
        {
            local_entry_count[i1] = 0;
        }
        else{

            if (read_offset + ceil((float)row_count[i1] / nprocs) > row_count[i1])
                local_entry_count[i1] = row_count[i1] - read_offset;
            else
                local_entry_count[i1] = (u64) ceil((float)row_count[i1] / nprocs);

            if (local_entry_count[i1] != 0)
            {

                char data_filename[1024];
                sprintf(data_filename, "%s/data.raw", argv[2 + 2*i1]);
                int fp = open(data_filename, O_RDONLY);


                size_t bsize = (u64) local_entry_count[i1] *  (u64)COL_COUNT * sizeof(u64);
                input_buffer[i1] = new u64[local_entry_count[i1] * COL_COUNT];


                u64 rb_size = pread(fp, input_buffer[i1], bsize, read_offset * COL_COUNT * sizeof(u64));
                if (rb_size != local_entry_count[i1] * COL_COUNT * sizeof(u64))
                {
                    std::cout << "Wrong IO: rank: " << rank << " " << rb_size << " " <<  bsize << " " << local_entry_count[i1] << std::endl;
                    MPI_Abort(MPI_COMM_WORLD, -1);
                }
                close(fp);
            }
        }

        ior_end[i1] = MPI_Wtime();
    }

    MPI_Barrier(MPI_COMM_WORLD);

    double start = MPI_Wtime();
    Relation1Map tempT;
    for (u32 i1 = 0; i1 < relation_count; i1++)
    {
        if (rank == 64)
            std::cout << "Local entry count () " << i1 << " " << local_entry_count[i1] << std::endl;

        hash_start[i1] = MPI_Wtime();
        /* process_size[j] stores the number of samples to be sent to process with rank j */
        int process_size[nprocs];
        memset(process_size, 0, nprocs * sizeof(int));

        /* vector[i] contains the data that needs to be sent to process i */
        std::vector<tuple<2>> *process_data_vector;
        process_data_vector = new std::vector<tuple<2>>[nprocs];

        tuple<2> dt;
        for (u32 i = 0; i < local_entry_count[i1] * COL_COUNT; i=i+2)
        {
            uint64_t index = outer_hash(input_buffer[i1][i])%nprocs;
            process_size[index] = process_size[index] + COL_COUNT;

            auto itx = tempT.find(input_buffer[i1][i]);
            if( itx != tempT.end() )
            {
                auto it2x = (itx->second)->find(input_buffer[i1][i + 1]);
                if( it2x != (itx->second)->end() )
                {
                }
                else
                {
                    (itx->second)->insert(std::make_pair(input_buffer[i1][i + 1], 0));
                    tempT[input_buffer[i1][i]] = itx->second;


                    dt[0] = input_buffer[i1][i];
                    dt[1] = input_buffer[i1][i + 1];

                    index = outer_hash(dt[0])%nprocs;
                    process_data_vector[index].push_back(dt);
                }
            }
            else
            {
                Relation0Map* k = new Relation0Map();
                k->insert(std::make_pair(input_buffer[i1][i + 1], 0));
                tempT[input_buffer[i1][i]] = k;

                index = outer_hash(input_buffer[i1][i])%nprocs;

                dt[0] = input_buffer[i1][i];
                dt[1] = input_buffer[i1][i + 1];
                process_data_vector[index].push_back(dt);
            }
        }

        hash_end[i1] = MPI_Wtime();

        comm_start[i1] = MPI_Wtime();

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
        hash_entry_count[i1] = recv_process_size_buffer[0];
        for (u32 i = 1; i < nprocs; i++)
        {
            prefix_sum_recv_process_size_buffer[i] = prefix_sum_recv_process_size_buffer[i - 1] + recv_process_size_buffer[i - 1];
            hash_entry_count[i1] = hash_entry_count[i1] + recv_process_size_buffer[i];
        }

        hashed_data[i1] = new u64[hash_entry_count[i1]];
        MPI_Alltoallv(process_data, process_size, prefix_sum_process_size, MPI_UNSIGNED_LONG_LONG, hashed_data[i1], recv_process_size_buffer, prefix_sum_recv_process_size_buffer, MPI_UNSIGNED_LONG_LONG, comm);

        hash_entry_count[i1] = hash_entry_count[i1] / COL_COUNT;

        delete[] process_data;

        comm_end[i1] = MPI_Wtime();


        union_start[i1] = MPI_Wtime();
        /*
        for ( u32 j1 = 0; j1 < hash_entry_count[i1] * COL_COUNT; j1 = j1 + 2)
        {
            total_tuple++;
            t1[0] = hashed_data[i1][j1];
            t1[1] = hashed_data[i1][j1+1];
            if (unionG.insert(t1))
                inserted_tuple++;
        }
        */

        for (u32 i = 0; i < hash_entry_count[i1] * COL_COUNT; i = i + 2)
        {
            total_tuple++;
            auto it = unionG.find(hashed_data[i1][i]);
            if ( it != unionG.end() ) {
                auto it2 = (it->second)->find(hashed_data[i1][i + 1]);
                if ( it2 != (it->second)->end() ) {
                }
                else {
                    (it->second)->insert(std::make_pair(hashed_data[i1][i + 1], 0));
                    unionG.insert(std::make_pair(hashed_data[i1][i], it->second));
                    inserted_tuple++;
                }
            }
            else {
                Relation0Map *k = new Relation0Map;
                k->insert(std::make_pair(hashed_data[i1][i + 1], 0));
                unionG.insert(std::make_pair(hashed_data[i1][i],k));
                inserted_tuple++;
            }
        }

        union_end[i1] = MPI_Wtime();
    }

    MPI_Barrier(MPI_COMM_WORLD);

    Relation1Map O;
    double join_start = MPI_Wtime();
    u64 jcount = parallel_join(O, unionG, unionG);
    double join_end = MPI_Wtime();

    double end = MPI_Wtime();

    double elapsed_time = end - start;
    double max_time = 0;
    u64 total_sum = 0;
    u64 total_sum2 = 0;
    MPI_Allreduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&inserted_tuple, &total_sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
    MPI_Allreduce(&total_tuple, &total_sum2, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);

    u64 join_total_sum = 0;
    MPI_Allreduce(&jcount, &join_total_sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);

    if (max_time == elapsed_time)
    {
        std::cout << nprocs << " : "
                  << "[L] Total tuple: " << total_tuple
                  << " [G] Total tuple: " << total_sum2
                  << " [L] Inserted tuple: " << inserted_tuple
                  << " [G] Inserted tuple: " << total_sum
                  << std::endl;

        double time = 0;

        double dedup = 0;
        double comm = 0;
        double insert = 0;
        for (u32 i = 0; i < relation_count; i++)
        {
            std::cout << "[" << i << "]: "
                      << " Global row count: " << row_count[i]
                      << " Local row count: " << local_entry_count[i]
                      << " Hash row count: " << hash_entry_count[i]
                      << " Read time: " << (ior_end[i] - ior_start[i])
                      << " Hash time: " << (hash_end[i] - hash_start[i])
                      << " comm time: " << (comm_end[i] - comm_start[i])
                      << " Insert time: " << (union_end[i] - union_start[i])
                      << std::endl;

            dedup = dedup + (hash_end[i] - hash_start[i]);
            comm = comm + (comm_end[i] - comm_start[i]);
            insert = insert + (union_end[i] - union_start[i]);

            time = time + (comm_end[i] - comm_start[i]) + (hash_end[i] - hash_start[i]) + (union_end[i] - union_start[i]);
        }

        std::cout << "Join Count " << join_total_sum << " " << join_end - join_start << std::endl;

        std::cout << "Read + hash + insert: " << time << std::endl;
        std::cout << "Total time: " << max_time << " Union " << (dedup + comm + insert) << " [ " << dedup << " " << comm << " " << insert << " ] Join [ " << (join_end - join_start)
                  << " " << (i2 - j1x) << " "
                  << (j2 - j1x) + (c2 - c1) + (i2 - in1) << " ] "
                  << (j2 - j1x)
                  << " " << (c2 - c1)
                  << " " << (i2 - in1)


                  << " " << " Total " << (dedup + comm + insert) + (join_end - join_start) << std::endl;
    }

    /*

    for (u32 i = 0; i < relation_count; i++)
    {
        delete[] input_buffer[i1];
        delete[] hashed_data[i1];
    }

    delete[] input_buffer;
    delete[] hashed_data;


    delete[]  row_count;
    delete[]  hash_entry_count;
    delete[]  local_entry_count;


    delete[] ior_start;
    delete[] ior_end;
    delete[] hash_start;
    delete[] hash_end;
    delete[] union_start;
    delete[] union_end;
    */

    MPI_Finalize();

    return 0;
}
