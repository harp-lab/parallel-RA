/*
 * copy
 * Copyright (c) Sidharth Kumar, et al, see License.md
 */



#include "../parallel_RA_inc.h"


void parallel_copy::local_copy(u32 buckets, google_relation** input, u32* input_bucket_map, relation* output, std::vector<int> reorder_map, u32 arity, u32 join_column_count, all_to_all_buffer& copy_buffer, int ra_counter)
{
    u32* output_sub_bucket_count = output->get_sub_bucket_per_bucket_count();
    u32** output_sub_bucket_rank = output->get_sub_bucket_rank();

    int projection_column_count=0;
    for (int x : reorder_map)
        if (x == -1)
            projection_column_count++;

    copy_buffer.width[ra_counter] = (output->get_arity() + 1) - projection_column_count;

    for (u32 i = 0; i < buckets; i++)
    {
        if (input_bucket_map[i] == 1)
        {
            input[i]->as_all_to_all_copy_buffer(copy_buffer, {}, reorder_map, ra_counter, buckets, output_sub_bucket_count, output_sub_bucket_rank, (arity+1), join_column_count);
        }
    }

    return;
}
