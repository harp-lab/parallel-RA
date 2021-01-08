/*
 * copy
 * Copyright (c) Sidharth Kumar, et al, see License.md
 */



#include "../parallel_RA_inc.h"

#ifdef GOOGLE_MAP
void parallel_copy_filter::local_copy_filter(u32 buckets, google_relation* input, u32* input_bucket_map, relation* output, std::vector<int> reorder_map, u32 arity, u32 join_column_count, all_to_allv_buffer& copy_filter_buffer, int ra_counter)
{
    u32* output_sub_bucket_count = output->get_sub_bucket_per_bucket_count();
    u32** output_sub_bucket_rank = output->get_sub_bucket_rank();


    copy_filter_buffer.width[ra_counter] = reorder_map.size();
    assert(copy_filter_buffer.width[ra_counter] == (int)output->get_arity());


    for (u32 i = 0; i < buckets; i++)
        if (input_bucket_map[i] == 1)
            input[i].as_all_to_allv_copy_filter_buffer(copy_filter_buffer, {}, reorder_map, ra_counter, buckets, output_sub_bucket_count, output_sub_bucket_rank, (arity), join_column_count, lambda, output->get_join_column_count(), output->get_is_canonical());



    return;
}


void parallel_copy_filter::local_copy_filter_with_threshold(int threshold, int ra_count,  u32 buckets, google_relation* input, u32* input_bucket_map, relation* output, std::vector<int> reorder_map, u32 arity, u32 join_column_count, all_to_all_buffer& copy_filter_buffer, int ra_counter)
{
    u32* output_sub_bucket_count = output->get_sub_bucket_per_bucket_count();
    u32** output_sub_bucket_rank = output->get_sub_bucket_rank();


    copy_filter_buffer.width[ra_counter] = reorder_map.size();
    assert(copy_filter_buffer.width[ra_counter] == (int)output->get_arity());


    for (u32 i = 0; i < buckets; i++)
        if (input_bucket_map[i] == 1)
            input[i].as_all_to_all_copy_filter_buffer(threshold, ra_count, copy_filter_buffer, {}, reorder_map, ra_counter, buckets, output_sub_bucket_count, output_sub_bucket_rank, (arity), join_column_count, lambda, output->get_join_column_count(), output->get_is_canonical());

    return;
}
#else


void parallel_copy_generate::local_copy_generate(u32 buckets, shmap_relation* input, u32* input_bucket_map, relation* output, std::vector<int> reorder_map, u32 arity, u32 join_column_count, all_to_allv_buffer& copy_filter_buffer, int ra_counter)
{
    u32* output_sub_bucket_count = output->get_sub_bucket_per_bucket_count();
    u32** output_sub_bucket_rank = output->get_sub_bucket_rank();


    copy_filter_buffer.width[ra_counter] = reorder_map.size();
    assert(copy_filter_buffer.width[ra_counter] == (int)output->get_arity());


    for (u32 i = 0; i < buckets; i++)
        if (input_bucket_map[i] == 1)
            input[i].as_all_to_allv_copy_generate_buffer(copy_filter_buffer, {}, reorder_map, ra_counter, buckets, output_sub_bucket_count, output_sub_bucket_rank, (arity), join_column_count, lambda, output->get_join_column_count(), output->get_is_canonical());

    return;
}
#endif
