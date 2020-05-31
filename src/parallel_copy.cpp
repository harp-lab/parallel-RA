#include "parallel_RA_inc.h"




void parallel_copy::local_copy(u32 buckets, google_relation* input, u32* input_bucket_map, relation* output, std::vector<int> reorder_map, vector_buffer* local_join_output, int* process_size, int* cumulative_process_size)
{
    u32 arity = output->get_arity();
    u32* output_sub_bucket_count = output->get_sub_bucket_per_bucket_count();
    u32** output_sub_bucket_rank = output->get_sub_bucket_rank();

    for (u32 i = 0; i < buckets; i++)
    {
        if (input_bucket_map[i] == 1)
        {
            vector_buffer temp_buffer;
            temp_buffer.vector_buffer_create_empty();

            std::vector<u64> prefix = {};
            input[i].as_vector_buffer_recursive(&temp_buffer, prefix);

            for (u32 s = 0; s < temp_buffer.size / sizeof(u64); s=s+arity)
            {
                u64 reordered_cur_path[arity];
                for (u32 j =0; j < arity; j++)
                    memcpy(reordered_cur_path + reorder_map[j], (&temp_buffer)->buffer + ((s + j) * sizeof(u64)), sizeof(u64));

                uint64_t bucket_id = hash_function(reordered_cur_path[0]) % buckets;
                uint64_t sub_bucket_id = hash_function(reordered_cur_path[1]) % output_sub_bucket_count[bucket_id];
                int index = output_sub_bucket_rank[bucket_id][sub_bucket_id];

                process_size[index] = process_size[index] + arity;
                cumulative_process_size[index] = cumulative_process_size[index] + arity;


                local_join_output[index].vector_buffer_append((const unsigned char*)reordered_cur_path, sizeof(u64)*arity);
            }
            temp_buffer.vector_buffer_free();
        }
    }

    return;
}

