// Class to manage basic MPI stuff

#ifndef COMM_H
#define COMM_H

#include "compat.h"
#include <mpi.h>

class mpi_comm
{

private:

    // stores the rank of a process
    int rank;

    // stores the total number of processes in a communicator
    int nprocs;


    // All processes within a comm will have the same number of buckets
    u32 buckets;

    // Stores the MPI comm
    MPI_Comm world_comm;

public:


    // sets the number of buckets = nprocs * bucket_factor
    void set_number_of_buckets(u32 bucket_factor)
    {
        buckets = bucket_factor * nprocs;
        std::cout << "Number of buckets " << buckets << std::endl;
    }


    // returns the total number of processes
    int get_nprocs()
    {
        return nprocs;
    }


    // returns the total number of buckets
    int get_number_of_buckets()
    {
        return buckets;
    }


    // returns the total number of processes
    int get_rank()
    {
        return rank;
    }


    // returns the associated communicator
    MPI_Comm get_comm()
    {
        return world_comm;
    }


    // initialize mpi
    void create(int argc, char **argv)
    {
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        world_comm = MPI_COMM_WORLD;

        std::cout << "[COMM] nprocs " << nprocs << " rank " << rank << std::endl;
    }


    // destroy mpi
    void destroy()
    {
        MPI_Finalize();
    }

};

#endif
