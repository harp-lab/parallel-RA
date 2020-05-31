#ifndef PARALLEL_COPY_H
#define PARALLEL_COPY_H



class parallel_copy: public parallel_RA
{

private:
    relation* copy_input0_table;
    int copy_input0_graph_type;
    relation* copy_output_table;

    std::vector<int> copy_reorder_index_array;

public:
    parallel_copy()
    {
        iteration_count = 1;
        RA_type = COPY;
    }

    parallel_copy(relation* G, int G_version, relation* output, std::vector<int> reorder_index_array)
        : copy_input0_table(G), copy_input0_graph_type(G_version), copy_output_table(output), copy_reorder_index_array(reorder_index_array)
    {
        RA_type = COPY;
    }



    void set_copy_input(relation* i0, int g_type0)
    {
        copy_input0_table = i0;
        copy_input0_graph_type = g_type0;
    }
    relation* get_copy_input()    {return copy_input0_table;}

    int get_copy_input0_graph_type() {return copy_input0_graph_type;}
    void set_copy_output(relation*& out) { copy_output_table = out; }
    relation* get_copy_output() {return copy_output_table; }

    void set_copy_rename_index (std::vector<int> pria) {copy_reorder_index_array = pria;}

    void get_copy_rename_index(std::vector<int>* projection_reorder_index_array) {*projection_reorder_index_array = this->copy_reorder_index_array;}


    void local_copy(u32 buckets, google_relation* input, u32* input_bucket_map, relation* output, std::vector<int> reorder_map, vector_buffer* local_join_output, int* process_size, int* cumulative_process_size);
};

#endif

