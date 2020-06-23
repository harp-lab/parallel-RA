#ifndef __all_to_all_H__
#define __all_to_all_H__


#if 1
struct all_to_all_buffer
{
    int ra_count;

    int nprocs;

    int *ra_size;

    // local_compute_output is of size RA_list size x nprocs, this contains data per RA rule for per process
    vector_buffer** local_compute_output;

    // local_compute_output_size is of size RA_list size x nprocs, this contains size of data per RA rule for per process
    int **local_compute_output_size;

    // cumulative_tuple_process_map is of size nprocs, cumulative_tuple_process_map[i] contains number of tuples to be transmitted to rank-i process, across all Relation Algebra (RA) rules
    int *cumulative_tuple_process_map;
};
#endif


void all_to_all_comm(vector_buffer* local_join_output, int* process_size, u64 *outer_hash_buffer_size, u64 **outer_hash_data, MPI_Comm comm);


#endif
