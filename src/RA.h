#ifndef RA_H
#define RA_H

#include <mpi.h>
#include "btree_relation.h"



class RA
{

public:

    relation<2> * parallel_join(relation<2>* delT, relation<2>& G, relation<2>& T, int lc, int* lb, int* running_t_count, double* running_time, u32 nprocs, u32 rank, MPI_Comm comm)
    {
        double j1, j2;
        double c1 = 0, c2 = 0;
        double i1 = 0, i2 = 0;
        double v1 = 0;
        double v2 = 0;

        j1 = MPI_Wtime();

        tuple<2> t;
        t[0] = -1;
        t[1] = -1;
        tuple<2> selectall(t);
        relation<2> tempT;

        // Send Join output
        /* process_size[j] stores the number of samples to be sent to process with rank j */
        int* process_size = new int[nprocs];
        memset(process_size, 0, nprocs * sizeof(int));

        /* vector[i] contains the data that needs to be sent to process i */
        std::vector<tuple<2>> *process_data_vector;
        process_data_vector = new std::vector<tuple<2>>[nprocs];

        u64 tuple_count = 0;


        tuple<2> dt;
        tuple<2> s;
        s[1] = -1;
        tuple<2> select;
        uint64_t index;
        for (relation<2>::iter dit(*delT, selectall); dit.more(); dit.advance())
        {
            s[0] = (*dit)[1];
            select = s;

            for (relation<2>::iter git(G, select); git.more(); git.advance())
            {
                dt[0] = (*dit)[0];
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


#if 1
        c1 = MPI_Wtime();
        int prefix_sum_process_size[nprocs];
        u64 non_deduplicate_tuple_count = process_data_vector[0].size();
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
        delete[] process_size;


        c2 = MPI_Wtime();

        i1 = MPI_Wtime();
        relation<2> * delTT = new relation<2>;
        u64 tcount = 0;
        u64 tduplicates = 0;
        for (u32 k = 0; k < outer_hash_buffer_size; k = k + 2)
        {
            tuple<2> dt;
            dt[0] = hash_buffer[k];
            dt[1] = hash_buffer[k + 1];

            if (T.insert(dt) == true)
            {
                tcount++;
                delTT->insert(dt);
            }
            else
                tduplicates++;
        }

        delete[] hash_buffer;
        delete[] process_data;
        i2 = MPI_Wtime();

        v1 = MPI_Wtime();
        if (lc % 10 == 0)
        {
            int sum = 0;
            MPI_Allreduce(&tcount, &sum, 1, MPI_INT, MPI_BOR, comm);
            if(sum == 0)
                *lb = 1;
            else
                *lb = 0;
        }
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

};

#endif



