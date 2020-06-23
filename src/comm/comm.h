// Class to manage basic MPI stuff

#ifndef COMM_H
#define COMM_H


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

    mpi_comm()
    {
    }

    mpi_comm (const mpi_comm &copy)
    {
        rank = copy.rank;
        nprocs = copy.nprocs;
        local_nprocs = copy.local_nprocs;
        local_rank = copy.local_rank;
        world_comm = copy.world_comm;
        local_comm = copy.local_comm;
    }



    int get_nprocs()    {return nprocs;}


    int get_local_nprocs()  {return local_nprocs;}


    int get_local_rank()    {return local_rank;}

    // returns the total number of processes
    int get_rank()  {return rank;}


    MPI_Comm get_comm() {return world_comm;}


    MPI_Comm get_local_comm()    {return local_comm;}


    void set_local_comm(MPI_Comm* comm);


    void create(int argc, char **argv);


    void destroy();


};

#endif
