/*
 * join
 * Copyright (c) Sidharth Kumar, et al, see License.md
 */



#include "../parallel_RA_inc.h"

void parallel_join::local_join(int join_order,
                               u32 buckets,
                               int input0_buffer_size, int input0_buffer_width, u64 *input0_buffer,
                               google_relation *input1, u32 i1_size, int input1_buffer_width,
                               std::vector<int> reorder_map_array,
                               relation* output,
                               all_to_all_buffer& join_buffer,
                               int counter,
                               int join_column_count,
                               u32* local_join_duplicates,
                               u32* local_join_inserts)
{
    join_buffer.width[counter] = reorder_map_array.size();//input0_buffer_width + input1_buffer_width - join_column_count - projection_column_count;


    google_relation deduplicate;
    u32* output_sub_bucket_count = output->get_sub_bucket_per_bucket_count();
    u32** output_sub_bucket_rank = output->get_sub_bucket_rank();



    if (join_order == LEFT)
    {
        for (int k1 = 0; k1 < input0_buffer_size; k1 = k1 + input0_buffer_width)
        {
            std::vector<u64> prefix;
            for (int jc=0; jc < join_column_count; jc++)
                prefix.push_back(input0_buffer[k1 + jc]);

            u64 bucket_id = tuple_hash(input0_buffer + k1, join_column_count) % buckets;

            //if (name == "test")
            //    std::cout << "LEFT Testing Join " << prefix[0] <<  std::endl;

            input1[bucket_id].as_all_to_all_left_join_buffer(prefix, join_buffer, input0_buffer + k1, input0_buffer_width, input1_buffer_width, counter, buckets, output_sub_bucket_count, output_sub_bucket_rank, reorder_map_array, join_column_count, deduplicate, local_join_duplicates, local_join_inserts);
        }
    }

    else if (join_order == RIGHT)
    {
        for (int k1 = 0; k1 < input0_buffer_size; k1 = k1 + input0_buffer_width)
        {
            std::vector<u64> prefix;
            for (int jc=0; jc < join_column_count; jc++)
                prefix.push_back(input0_buffer[k1 + jc]);

            u64 bucket_id = tuple_hash(input0_buffer + k1, join_column_count) % buckets;

#if 0
            if (name == "test")
                std::cout << "RIGHT Testing Join Rel 1 length "
                          << input0_buffer_width << ", Rel 2 length "
                          << input1_buffer_width << ", VAL "
                          << input0_buffer[k1] << ", "
                          << input0_buffer[k1 + 1] << ", "
                          << input0_buffer[k1 + 2] << ", "
                          << input0_buffer[k1 + 3]
                          <<  std::endl;
#endif

            input1[bucket_id].as_all_to_all_right_join_buffer(prefix, join_buffer, input0_buffer + k1, input0_buffer_width, input1_buffer_width, counter, buckets, output_sub_bucket_count, output_sub_bucket_rank, reorder_map_array, join_column_count, deduplicate, local_join_duplicates, local_join_inserts, "test");
        }
    }

    deduplicate.remove_tuple();
    return;
}
