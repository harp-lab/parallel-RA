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


    // Stores the MPI comm
    MPI_Comm world_comm;

public:

    // returns the total number of processes
    int get_nprocs()
    {
        return nprocs;
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

        //std::cout << "[COMM] nprocs " << nprocs << " rank " << rank << std::endl;
    }


    // destroy mpi
    void destroy()
    {
        MPI_Finalize();
    }

};

#endif
