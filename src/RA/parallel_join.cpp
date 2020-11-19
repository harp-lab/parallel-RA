/*
 * join
 * Copyright (c) Sidharth Kumar, et al, see License.md
 */



#include "../parallel_RA_inc.h"

bool parallel_join::local_join(int threshold, int* offset,
                               int join_order,
                               u32 buckets,
                               int input0_buffer_size, int input0_buffer_width, u64 *input0_buffer,
                               google_relation *input1, u32 i1_size, int input1_buffer_width,
                               std::vector<int> reorder_map_array,
                               relation* output,
                               all_to_allv_buffer& join_buffer,
                               int counter,
                               int join_column_count,
                               u32* global_join_duplicates,
                               u32* global_join_inserts)
{
    join_buffer.width[counter] = reorder_map_array.size();

    google_relation deduplicate;
    u32* output_sub_bucket_count = output->get_sub_bucket_per_bucket_count();
    u32** output_sub_bucket_rank = output->get_sub_bucket_rank();

    if (*offset > input0_buffer_size || input0_buffer_size == 0 || i1_size == 0)
        return true;

    int local_join_count=0;
    if (join_order == LEFT)
    {
        for (int k1 = *offset; k1 < input0_buffer_size; k1 = k1 + input0_buffer_width)
        {
            std::vector<u64> prefix;
            for (int jc=0; jc < join_column_count; jc++)
                prefix.push_back(input0_buffer[k1 + jc]);

            u64 bucket_id = tuple_hash(input0_buffer + k1, join_column_count) % buckets;

            input1[bucket_id].as_all_to_allv_left_join_buffer(prefix, join_buffer, input0_buffer + k1, input0_buffer_width, input1_buffer_width, counter, buckets, output_sub_bucket_count, output_sub_bucket_rank, reorder_map_array, join_column_count, deduplicate, &local_join_count, global_join_duplicates, global_join_inserts, output->get_join_column_count(), output->get_is_canonical());

            //std::cout << "local_join_count " << local_join_count << " Threshold " << threshold << " k1 " << k1 << " offset " << *offset << " " << input0_buffer_width << std::endl;
            if (local_join_count > threshold)
            {
                *offset = k1 + input0_buffer_width;
                deduplicate.remove_tuple();
                return false;
            }
        }
    }

    else if (join_order == RIGHT)
    {
        for (int k1 = *offset; k1 < input0_buffer_size; k1 = k1 + input0_buffer_width)
        {
            std::vector<u64> prefix;
            for (int jc=0; jc < join_column_count; jc++)
                prefix.push_back(input0_buffer[k1 + jc]);

            u64 bucket_id = tuple_hash(input0_buffer + k1, join_column_count) % buckets;

            input1[bucket_id].as_all_to_allv_right_join_buffer(prefix, join_buffer, input0_buffer + k1, input0_buffer_width, input1_buffer_width, counter, buckets, output_sub_bucket_count, output_sub_bucket_rank, reorder_map_array, join_column_count, deduplicate, &local_join_count, global_join_duplicates, global_join_inserts, "test", output->get_join_column_count(), output->get_is_canonical());

            //std::cout << "local_join_count " << local_join_count << " Threshold " << threshold << " k1 " << k1 << " offset " << *offset << " " << input0_buffer_width << std::endl;
            if (local_join_count > threshold)
            {
                *offset = k1 + input0_buffer_width;
                //std::cout << "Setting offset " << *offset << std::endl;
                deduplicate.remove_tuple();
                return false;
            }
        }
    }

    deduplicate.remove_tuple();
    return true;
}


bool parallel_join::local_join_with_threshold(int threshold, int RA_count, int* offset,
                               int join_order,
                               u32 buckets,
                               int input0_buffer_size, int input0_buffer_width, u64 *input0_buffer,
                               google_relation *input1, u32 i1_size, int input1_buffer_width,
                               std::vector<int> reorder_map_array,
                               relation* output,
                               all_to_all_buffer& join_buffer,
                               int counter,
                               int join_column_count,
                               u32* global_join_duplicates,
                               u32* global_join_inserts)
{
    join_buffer.width[counter] = reorder_map_array.size();

    google_relation deduplicate;
    u32* output_sub_bucket_count = output->get_sub_bucket_per_bucket_count();
    u32** output_sub_bucket_rank = output->get_sub_bucket_rank();

    if (*offset > input0_buffer_size || input0_buffer_size == 0 || i1_size == 0)
        return true;

    u64 *temp_buffer = new u64[i1_size * reorder_map_array.size()];
    int local_join_count=0;
    if (join_order == LEFT)
    {
        for (int k1 = *offset; k1 < input0_buffer_size; k1 = k1 + input0_buffer_width)
        {
            std::vector<u64> prefix;
            for (int jc=0; jc < join_column_count; jc++)
                prefix.push_back(input0_buffer[k1 + jc]);

            u64 bucket_id = tuple_hash(input0_buffer + k1, join_column_count) % buckets;

            u32 tracker = 0;
            input1[bucket_id].as_all_to_all_left_join_buffer(threshold, RA_count, prefix, join_buffer, input0_buffer + k1, input0_buffer_width, input1_buffer_width, counter, buckets, output_sub_bucket_count, output_sub_bucket_rank, reorder_map_array, join_column_count, deduplicate, &local_join_count, global_join_duplicates, global_join_inserts, output->get_join_column_count(), output->get_is_canonical(), &tracker);
            if (tracker == 1)
            {
                *offset = k1;
                deduplicate.remove_tuple();
                return false;
            }
        }
    }

    else if (join_order == RIGHT)
    {
        for (int k1 = *offset; k1 < input0_buffer_size; k1 = k1 + input0_buffer_width)
        {
            std::vector<u64> prefix;
            for (int jc=0; jc < join_column_count; jc++)
                prefix.push_back(input0_buffer[k1 + jc]);

            u64 bucket_id = tuple_hash(input0_buffer + k1, join_column_count) % buckets;

            u32 tracker = 0;
            input1[bucket_id].as_all_to_all_right_join_buffer(threshold, RA_count, prefix, join_buffer, input0_buffer + k1, input0_buffer_width, input1_buffer_width, counter, buckets, output_sub_bucket_count, output_sub_bucket_rank, reorder_map_array, join_column_count, deduplicate, &local_join_count, global_join_duplicates, global_join_inserts, "test", output->get_join_column_count(), output->get_is_canonical(), &tracker);
            if (tracker == 1)
            {
                *offset = k1;
                deduplicate.remove_tuple();
                return false;
            }
        }
    }

    delete[] temp_buffer;
    deduplicate.remove_tuple();
    return true;
}
