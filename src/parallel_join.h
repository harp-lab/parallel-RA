#ifndef PARALLEL_JOIN_H
#define PARALLEL_JOIN_H



class parallel_join: public parallel_RA
{

private:

    relation* join_input0_table;
    int join_input0_graph_type;

    relation* join_input1_table;
    int join_input1_graph_type;

    relation* join_output_table;

    int join_column_count;

    std::vector<int> projection_reorder_index_array;
    int projection_reorder_index_array_length;

public:
    parallel_join()
    {
        RA_type = JOIN;
    }

    parallel_join(relation* G, int G_type, relation* T, int T_type, relation* output, int jc_count, std::vector<int> projection_reorder_index_array)
        : join_input0_table(G), join_input0_graph_type(G_type), join_input1_table(T), join_input1_graph_type(T_type), join_output_table(output), join_column_count(jc_count), projection_reorder_index_array(projection_reorder_index_array)
    {
        RA_type = JOIN;
    }


    void set_join_input0(relation* i0, int g_type0)
    {
        join_input0_table = i0;
        join_input0_graph_type = g_type0;
    }

    relation* get_join_input0() {return join_input0_table;}

    int get_join_input0_graph_type()    {return join_input0_graph_type;}

    void set_join_input1(relation* i1, int g_type1)
    {
        join_input1_table = i1;
        join_input1_graph_type = g_type1;
    }

    relation* get_join_input1() {return join_input1_table;}

    int get_join_input1_graph_type()    {return join_input1_graph_type;}

    void set_join_output(relation*& out)    {join_output_table = out;}

    relation* get_join_output() {return join_output_table;}

    void set_join_projection_index (std::vector<int> pria)  {projection_reorder_index_array = pria;}

    void set_join_column_count (int jcc)    {join_column_count = jcc;}

    int get_join_column_count ()    {return join_column_count;}

    void get_join_projection_index(std::vector<int>* projection_reorder_index_array)    {*projection_reorder_index_array = this->projection_reorder_index_array; }

    // TODO join_order
    u32 local_join( u32 buckets,
                    int input0_buffer_size, int input0_buffer_width, u64 *input0_buffer, int join_order,
                    google_relation *input1, u32 i1_size, int input1_buffer_width,
                    std::vector<int> reorder_map_array,
                    relation* output,
                    vector_buffer** local_join_output, int** process_size, int* cumulative_process_size,
                    u32 threshhold, int join_colun_count,
                    u32* local_join_count, int iteration);

    void intra_bucket_comm(u32 buckets,
                           google_relation *rel,
                           int* input_distinct_sub_bucket_rank_count, int** input_distinct_sub_bucket_rank, u32* input_bucket_map,
                           int* output_distinct_sub_bucket_rank_count, int** output_distinct_sub_bucket_rank, u32* output_bucket_map,
                           u64 *total_buffer_size, u64 **recvbuf);

};


#endif
