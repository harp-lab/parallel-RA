#ifndef INTRA_BUCKET_COMM_H
#define INTRA_BUCKET_COMM_H

void intra_bucket_comm(u32 buckets,
                       google_relation *rel,
                       int* input_distinct_sub_bucket_rank_count, int** input_distinct_sub_bucket_rank, u32* input_bucket_map,
                       int* output_distinct_sub_bucket_rank_count, int** output_distinct_sub_bucket_rank, u32* output_bucket_map,
                       u64 *total_buffer_size, u64 **recvbuf,
                       MPI_Comm mcomm);

#endif
