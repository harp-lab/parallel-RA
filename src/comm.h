// Class to manage basic MPI stuff

#ifndef COMM_H
#define COMM_H

#include "compat.h"
#include <mpi.h>

class mpi_comm
{

private:

    //u32 buckets;

    // stores the rank of a process
    int rank;

    // stores the total number of processes in a communicator
    int nprocs;

    int local_nprocs;

    int local_rank;

    // Stores the MPI comm
    MPI_Comm world_comm;
    MPI_Comm local_comm;

public:

    mpi_comm ()
    {

    }

    mpi_comm (const mpi_comm &copy)
    {
        //buckets = copy.buckets;
        rank = copy.rank;
        nprocs = copy.nprocs;
        local_nprocs = copy.local_nprocs;
        local_rank = copy.local_rank;
        world_comm = copy.world_comm;
        local_comm = copy.local_comm;
    }

    // returns the total number of processes
    int get_nprocs()
    {
        return nprocs;
    }

    //int get_bucket_count()
    //{
    //    return buckets;
    //}


    //void get_bucket_count(int buck)
    //{
    //    buckets = buck;
    //}


    int get_local_nprocs()
    {

        return local_nprocs;
    }


    int get_local_rank()
    {

        return local_rank;
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


    MPI_Comm get_local_comm()
    {
        return local_comm;
    }


    void set_local_comm(MPI_Comm* comm)
    {
        local_comm = *comm;
        MPI_Comm_size(local_comm, &local_nprocs);
        MPI_Comm_rank(local_comm, &local_rank);
        //buckets = local_nprocs;
    }


    // initialize mpi
    void create(int argc, char **argv)
    {
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        world_comm = MPI_COMM_WORLD;

        //std::cout << "[COMM] nprocs " << nprocs << " rank " << rank << std::endl;
    }


    // destroy mpi
    void destroy()
    {
        MPI_Finalize();
    }

};

#endif
