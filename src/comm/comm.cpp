#include "../parallel_RA_inc.h"





void mpi_comm::set_local_comm(MPI_Comm* comm)
{
    local_comm = *comm;
    MPI_Comm_size(local_comm, &local_nprocs);
    MPI_Comm_rank(local_comm, &local_rank);
}



void mpi_comm::create(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    world_comm = MPI_COMM_WORLD;
}



void mpi_comm::destroy()
{
    MPI_Finalize();
}
