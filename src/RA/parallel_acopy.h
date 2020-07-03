/*
 * acopy
 * Copyright (c) Sidharth Kumar, et al, see License.md
 */



#ifndef PARALLEL_ACOPY_H
#define PARALLEL_ACOPY_H



class parallel_acopy: public parallel_RA
{

private:

    relation* acopy_input0_table;
    int acopy_input0_graph_type;

    relation* acopy_output_table;

    std::vector<int> acopy_reorder_index_array;

public:
    parallel_acopy()
    {
        RA_type = ACOPY;
    }

    parallel_acopy(relation* dest, relation* src, int src_version, std::vector<int> reorder_index_array)
        : acopy_input0_table(src), acopy_input0_graph_type(src_version), acopy_output_table(dest), acopy_reorder_index_array(reorder_index_array)
    {
        RA_type = ACOPY;
    }

    relation* get_acopy_input(){return acopy_input0_table;}
    relation* get_acopy_output() {return acopy_output_table;}
    int get_acopy_input0_graph_type() {return acopy_input0_graph_type;}
    void get_acopy_rename_index(std::vector<int>* projection_reorder_index_array) {*projection_reorder_index_array = this->acopy_reorder_index_array;}
    void local_acopy(u32 buckets, google_relation** input, u32* input_bucket_map, relation* output, std::vector<int> reorder_map, u32 arity, u32 join_column_count, all_to_all_buffer& acopy_buffer, int ra_counter);
};

#endif

