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

#define COL_COUNT 2

#include "btree/btree_map.h"

char TDname[1024];
typedef btree::btree_map<u64, u64> Relation0Map;
typedef btree::btree_map<u64, btree::btree_map<u64, u64>* > Relation1Map;

static int rank = 0;
static int nprocs = 1;
static MPI_Comm comm;
static u32 global_row_count;



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

void buffer_data_to_hash_buffer(u32 local_number_of_rows, u64* input_data,  int hash_column_index, u64** outer_hash_data, u32* outer_hash_buffer_size, MPI_Comm comm, int buckets, u32* subbuckets_G)
{
    /* process_size[j] stores the number of samples to be sent to process with rank j */
    int process_size[nprocs];
    memset(process_size, 0, nprocs * sizeof(int));

    //int tprocess_size[nprocs];
    //memset(tprocess_size, 0, nprocs * sizeof(int));

    /* vector[i] contains the data that needs to be sent to process i */
    std::vector<u64> process_data_vector[nprocs];

    if (hash_column_index == 0)
    {
        for (u32 i = 0; i < local_number_of_rows * COL_COUNT; i=i+2)
        {
            uint64_t bucket_id = outer_hash(input_data[i]) % buckets;
            uint64_t sub_bucket_id = outer_hash(input_data[i+1]) % subbuckets_G[bucket_id];

            //tprocess_size[bucket_id] = tprocess_size[bucket_id] + COL_COUNT;

            uint64_t index = (outer_hash((bucket_id << 32)^sub_bucket_id))%nprocs;
            process_size[index] = process_size[index] + COL_COUNT;

            process_data_vector[index].push_back(input_data[i]);
            process_data_vector[index].push_back(input_data[i + 1]);
        }
    }

    int prefix_sum_process_size[nprocs];
    memset(prefix_sum_process_size, 0, nprocs * sizeof(int));


    //if (rank == 0)
    //std::cout << "Dist " << 0 << " " << tprocess_size[0] << " " << process_size[0] << std::endl;
    for (int i = 1; i < nprocs; i++)
    {
        prefix_sum_process_size[i] = prefix_sum_process_size[i - 1] + process_size[i - 1];
        //if (rank == 0)
        //    std::cout << "Dist " << i << " " << tprocess_size[i] << " " << process_size[i] << std::endl;
    }

    /*
    if (rank == 0)
        for (int i = 0; i < nprocs; i++)
            std::cout << "i " << i << " : " << process_size[i]  << std::endl;
    */

    int process_data_buffer_size = prefix_sum_process_size[nprocs - 1] + process_size[nprocs - 1];

    //u64 process_data[process_data_buffer_size];
    u64* process_data = new u64[process_data_buffer_size];

    for (int i = 0; i < nprocs; i++)
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
    for (int i = 1; i < nprocs; i++)
        prefix_sum_recv_process_size_buffer[i] = prefix_sum_recv_process_size_buffer[i - 1] + recv_process_size_buffer[i - 1];


    /* Sending data to all processes */
    /* What is the buffer size to allocate */
    *outer_hash_buffer_size = 0;
    for(int i = 0; i < nprocs; i++)
        *outer_hash_buffer_size = *outer_hash_buffer_size + recv_process_size_buffer[i];



    uint total_row_size = 0;
    MPI_Allreduce(outer_hash_buffer_size, &total_row_size, 1, MPI_INT, MPI_SUM, comm);
    if(total_row_size != global_row_count * COL_COUNT)
    {
        printf("Incorrect distribution %d != %d %d\n", total_row_size, global_row_count, COL_COUNT);
        MPI_Abort(comm, -1);
    }
#if 1
    *outer_hash_data = new u64[*outer_hash_buffer_size];
    MPI_Alltoallv(process_data, process_size, prefix_sum_process_size, MPI_UNSIGNED_LONG_LONG, *outer_hash_data, recv_process_size_buffer, prefix_sum_recv_process_size_buffer, MPI_UNSIGNED_LONG_LONG, comm);

    *outer_hash_buffer_size = *outer_hash_buffer_size / COL_COUNT;

    delete[] process_data;
#endif
    return;
}




void load_balance_G(u32 buckets, u32** gmap_sub_bucket, u32* subbuckets_G, Relation1Map **G, std::ofstream& myfile)
{
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

    //if (rank == 1)
    //    for (u32 i = 0; i < buckets; i++)
    //        for (u32 j = 0; j < subbuckets_G[i]; j++)
    //            std::cout << "[Before] Element count at " << i << " " << j << " is " << gmap_sub_bucket[i][j] << std::endl;

    //if (rank == 0)
    myfile << "Rank " << rank << " max_sub_bucket_size " << max_sub_bucket_size << " min_sub_bucket_size " << min_sub_bucket_size << "\n";

    int global_max;
    int global_min;
    u64 total_outer_hash_buffer_size = 0;
    MPI_Allreduce(&max_sub_bucket_size, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&min_sub_bucket_size, &global_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);


    int* g_new_sub_bucket = new int[buckets];
    memset(g_new_sub_bucket, 0, buckets * sizeof(int));

    int* global_g_new_sub_bucket = new int[buckets];
    memset(global_g_new_sub_bucket, 0, buckets * sizeof(int));

    for (u32 i = 0; i < buckets; i++)
        for (u32 j = 0; j < subbuckets_G[i]; j++)
            g_new_sub_bucket[i] = g_new_sub_bucket[i] + gmap_sub_bucket[i][j] / global_min;


    MPI_Allreduce(g_new_sub_bucket, global_g_new_sub_bucket, buckets, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    //if (rank == 1)
    //    for (u32 b = 0; b < buckets; b++)
    //        std::cout << "b " << b << " " << subbuckets_G[b] << " " << global_g_new_sub_bucket[b] << std::endl;

    for (u32 i = 0; i < buckets; i++)
        delete[] gmap_sub_bucket[i];

    for (u32 ix = 0; ix < buckets; ix++)
    {
        gmap_sub_bucket[ix] = new u32[global_g_new_sub_bucket[ix]];
        memset(gmap_sub_bucket[ix], 0, sizeof(u32) * global_g_new_sub_bucket[ix]);
    }


    for (u32 i = 0; i < buckets; i++)
    {
#if 1
        /* process_size[j] stores the number of samples to be sent to process with rank j */
        int process_size[nprocs];
        memset(process_size, 0, nprocs * sizeof(int));

        /* vector[i] contains the data that needs to be sent to process i */
        std::vector<u64> process_data_vector[nprocs];

        int prefix_sum_process_size[nprocs];
        memset(prefix_sum_process_size, 0, nprocs * sizeof(int));

        int recv_process_size_buffer[nprocs];
        memset(recv_process_size_buffer, 0, nprocs * sizeof(int));

        //if (global_g_new_sub_bucket[i] != subbuckets_G[i])
        //{
        for (u32 j = 0; j < subbuckets_G[i]; j++)
        {
            for (auto it = G[i][j].begin(); it != G[i][j].end(); it++)
            {
                Relation0Map* it2 = it->second;
                for (auto dit2 = it2->begin(); dit2 != it2->end(); dit2++)
                {
                    uint64_t bucket_id = outer_hash(it->first) % buckets;
                    uint64_t sub_bucket_id = outer_hash(dit2->first) % global_g_new_sub_bucket[bucket_id];

                    uint64_t index = (outer_hash((bucket_id << 32)^sub_bucket_id))%nprocs;
                    process_size[index] = process_size[index] + COL_COUNT;

                    process_data_vector[index].push_back(it->first);
                    process_data_vector[index].push_back(dit2->first);
                }
            }
        }


        //if (rank == 1)
        //    std::cout << "process_size[0] " << process_size[0] << " process_size[1] " << process_size[1] << std::endl;

        for (int ix = 1; ix < nprocs; ix++)
            prefix_sum_process_size[ix] = prefix_sum_process_size[ix - 1] + process_size[ix - 1];

        int process_data_buffer_size = prefix_sum_process_size[nprocs - 1] + process_size[nprocs - 1];

        //u64 process_data[process_data_buffer_size];
        u64* process_data = new u64[process_data_buffer_size];

        for (int ix = 0; ix < nprocs; ix++)
            memcpy(process_data + prefix_sum_process_size[ix], &process_data_vector[ix][0], process_data_vector[ix].size() * sizeof(u64));


        /* This step prepares for actual data transfer */
        /* Every process sends to every other process the amount of data it is going to send */


        //if (rank == 1)
        //    for (u32 x = 0; x < nprocs; x++)
        //        std::cout << "B " << i << " process_size" << process_size[x] << std::endl;
        //}

        MPI_Alltoall(process_size, 1, MPI_INT, recv_process_size_buffer, 1, MPI_INT, MPI_COMM_WORLD);


        int prefix_sum_recv_process_size_buffer[nprocs];
        memset(prefix_sum_recv_process_size_buffer, 0, nprocs * sizeof(int));
        for (int ix = 1; ix < nprocs; ix++)
            prefix_sum_recv_process_size_buffer[ix] = prefix_sum_recv_process_size_buffer[ix - 1] + recv_process_size_buffer[ix - 1];


        /* Sending data to all processes */
        /* What is the buffer size to allocate */
        int outer_hash_buffer_size = 0;
        for(int ix = 0; ix < nprocs; ix++)
            outer_hash_buffer_size = outer_hash_buffer_size + recv_process_size_buffer[ix];


        //std::cout << "Rank " << rank << " bucket " << i << " size " << outer_hash_buffer_size << std::endl;


        u64* outer_hash_data = new u64[outer_hash_buffer_size];
        MPI_Alltoallv(process_data, process_size, prefix_sum_process_size, MPI_UNSIGNED_LONG_LONG, outer_hash_data, recv_process_size_buffer, prefix_sum_recv_process_size_buffer, MPI_UNSIGNED_LONG_LONG, comm);

        delete[] process_data;


        for (u32 j = 0; j < subbuckets_G[i]; j++)
        {
            Relation1Map::iterator iy2 = G[i][j].begin();
            for(; iy2 != G[i][j].end(); iy2++)
                delete (iy2->second);

            G[i][j].clear();
        }

        G[i] = new Relation1Map[global_g_new_sub_bucket[i]];

        total_outer_hash_buffer_size = total_outer_hash_buffer_size + outer_hash_buffer_size/2;

        //if (rank == 1)
        //std::cout << "outer_hash_buffer_size " << outer_hash_buffer_size << std::endl;
        //subbuckets_G[i] = global_g_new_sub_bucket[i];
        for (int in = 0; in < outer_hash_buffer_size; in = in + 2)
        {
            uint64_t bucket_id = outer_hash((outer_hash_data)[in]) % buckets;
            uint64_t sub_bucket_id = outer_hash((outer_hash_data)[in+1]) % global_g_new_sub_bucket[bucket_id];
            gmap_sub_bucket[bucket_id][sub_bucket_id]++;

            //if (rank == 1)
            //std::cout << "Val " << outer_hash_data[in] << " b " << bucket_id << " sbi " << sub_bucket_id << std::endl;

            auto it = G[bucket_id][sub_bucket_id].find((outer_hash_data)[in]);
            if( it != G[bucket_id][sub_bucket_id].end() ) {
                auto it2 = (it->second)->find((outer_hash_data)[in+1]);
                if( it2 != (it->second)->end() ) {
                    ;
                }
                else{
                    (it->second)->insert(std::make_pair((outer_hash_data)[in + 1], 0));
                    G[bucket_id][sub_bucket_id][(outer_hash_data)[in]] = it->second;
                }
            }
            else {
                Relation0Map *k = new Relation0Map;
                k->insert(std::make_pair((outer_hash_data)[in + 1], 0));
                G[bucket_id][sub_bucket_id].insert(std::make_pair((outer_hash_data)[in],k));
            }
        }

        delete[] outer_hash_data;
#endif

    }

#if 1
    //
    if (rank == 0)
    {
        for (u32 i = 0; i < buckets; i++)
            myfile << i << " Bucket count changed from " << subbuckets_G[i] << " to " << global_g_new_sub_bucket[i] << "\n";

        //for (u32 i = 0; i < buckets; i++)
        //    for (u32 j = 0; j < global_g_new_sub_bucket[i]; j++)
        //        std::cout << "[After] Element count at " << i << " " << j << " is " << gmap_sub_bucket[i][j] << std::endl;
    }
    //

    memcpy(subbuckets_G, global_g_new_sub_bucket, sizeof(u32)*buckets);

    int send_buf = total_outer_hash_buffer_size;
    int recv_buf_max = 0;
    int recv_buf_min = 0;
    MPI_Allreduce(&send_buf, &recv_buf_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&send_buf, &recv_buf_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    if (rank == 0)
    std::cout << "[Balanced]    Buckets " << buckets
              << " Sub-buckets " << subbuckets_G[0]
              << " Maximum buffer size " << recv_buf_max
              << " Minimum buffer size " << recv_buf_min
              << " Ratio " << (float)recv_buf_max/recv_buf_min
              << std::endl
              << std::endl;

    myfile << "Balanced [G] Rank: " << rank << " Hashed count: " << total_outer_hash_buffer_size << " Ratio: " << (float)((float) total_outer_hash_buffer_size/global_row_count)*100 << "\n";
#endif
}



void one_iteration(u32 local_row_count, u32 buckets, u32 *subbuckets_G, u64 *input_buffer)
{
    u32 G_hash_entry_count;
    u64 *G_hashed_data = NULL;
    buffer_data_to_hash_buffer(local_row_count, input_buffer, 0, &G_hashed_data, &G_hash_entry_count, MPI_COMM_WORLD, buckets, subbuckets_G);

    std::ofstream myfile;
    char filename[1024];
    sprintf(filename, "%s/filename_%d_%d_%d_%d", TDname, nprocs, buckets, subbuckets_G[0], rank);
    myfile.open (filename);


    myfile << "[G] Rank: " << rank << "Buckets: " << buckets << " Subbuckets: " << subbuckets_G[0] << " Hashed count: " << G_hash_entry_count << " Ratio: " << (float)((float) G_hash_entry_count/global_row_count)*100 << "\n";

    int send_buf = G_hash_entry_count;
    int recv_buf_max = 0;
    int recv_buf_min = 0;

    MPI_Allreduce(&send_buf, &recv_buf_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&send_buf, &recv_buf_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    if (rank == 0)
    std::cout << "[Un-balanced] Buckets " << buckets
              << " Sub-buckets " << subbuckets_G[0]
              << " Maximum buffer size " << recv_buf_max
              << " Minimum buffer size " << recv_buf_min
              << " Ratio " << (float)recv_buf_max/recv_buf_min
              << std::endl;


    u32 **gmap_sub_bucket = new u32*[buckets];
    memset(gmap_sub_bucket, 0, sizeof(u32*) * buckets);

    for (u32 i = 0; i < buckets; i++)
    {
        gmap_sub_bucket[i] = new u32[subbuckets_G[i]];
        memset(gmap_sub_bucket[i], 0, sizeof(u32) * subbuckets_G[i]);
    }

    Relation1Map **G = new Relation1Map*[buckets];
    for (u32 i = 0; i < buckets; i++)
        G[i] = new Relation1Map[subbuckets_G[i]];

    for (u32 i = 0; i < COL_COUNT * G_hash_entry_count; i = i + 2)
    {
        uint64_t bucket_id = outer_hash(G_hashed_data[i]) % buckets;
        uint64_t sub_bucket_id = outer_hash(G_hashed_data[i+1]) % subbuckets_G[bucket_id];
        gmap_sub_bucket[bucket_id][sub_bucket_id]++;

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

    load_balance_G(buckets, gmap_sub_bucket, subbuckets_G, G, myfile);

    for (u32 i = 0; i < buckets; i++)
    {
        for (u32 j = 0; j < subbuckets_G[i]; j++)
        {
            Relation1Map::iterator iy2 = G[i][j].begin();
            for(; iy2 != G[i][j].end(); iy2++)
                delete (iy2->second);
        }

        delete[] G[i];

    }

    myfile.close();

    delete[] G;
    delete[] gmap_sub_bucket;
}



//create big 2d array bucket/sub-bucket most will be empty
int main(int argc, char **argv)
{
    // Initializing MPI
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    comm = MPI_COMM_WORLD;

    u32 buckets;

    // this is the initialization step, we can change this as well
    buckets = nprocs * 1;//ceil((float)(3 * nprocs) / 4);

    // the relations can have different number of sub buckets


    u32 entry_count;
    u64 *input_buffer = NULL;
    parallel_read_input_relation_from_file_to_local_buffer(argv[1], &input_buffer, &entry_count);



    sprintf(TDname, "%d_%s", nprocs, argv[1]);
    if (rank == 0)
        mkdir(TDname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    MPI_Barrier(MPI_COMM_WORLD);

    for (u32 b = 1; b < 6; b++)
    {
        buckets = nprocs * b;

        u32 *subbuckets_G = new u32[buckets];

        for (u64 b = 0; b < buckets; b++)
            subbuckets_G[b] = 1;
        one_iteration(entry_count, buckets, subbuckets_G, input_buffer);

        for (u64 b = 0; b < buckets; b++)
            subbuckets_G[b] = 5;
        one_iteration(entry_count, buckets, subbuckets_G, input_buffer);

        for (u64 b = 0; b < buckets; b++)
            subbuckets_G[b] = 10;
        one_iteration(entry_count, buckets, subbuckets_G, input_buffer);

        for (u64 b = 0; b < buckets; b++)
            subbuckets_G[b] = 15;
        one_iteration(entry_count, buckets, subbuckets_G, input_buffer);

        for (u64 b = 0; b < buckets; b++)
            subbuckets_G[b] = 20;
        one_iteration(entry_count, buckets, subbuckets_G, input_buffer);

        delete[] subbuckets_G;
    }




    delete[] input_buffer;
    // Finalizing MPI
    MPI_Finalize();

    return 0;
}
