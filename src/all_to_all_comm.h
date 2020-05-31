#ifndef __all_to_all_H__
#define __all_to_all_H__


void all_to_all_comm(vector_buffer* local_join_output, int* process_size, u64 *outer_hash_buffer_size, u64 **outer_hash_data, MPI_Comm comm);


#endif
