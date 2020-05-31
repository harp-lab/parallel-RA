#ifndef RAM_H
#define RAM_H



class RAM
{
private:



    bool enable_dump_io = false;
    int iteration_count = -1;

    bool logging = false;
    bool comm_compaction=false;
    u32 threshold = 1000000;

    int refinement_chooser = 0;
    double refinement_factor = 4;

    u32 refinement_ts=10;

    std::vector<relation*> relation_manager;
    std::vector<parallel_RA*> RA_list;

    u64 *intra_bucket_buf_output_size;
    u64 **intra_bucket_buf_output;

    // local_compute_output is of size RA_list size x nprocs
    // this contains data per RA rule for per process
    vector_buffer** local_compute_output;

    // local_compute_output_size is of size RA_list size x nprocs
    // this contains size of data per RA rule for per process
    int **local_compute_output_size;

    // cumulative_tuple_process_map is of size nprocs
    // cumulative_tuple_process_map[i] contains number of tuples to be transmitted to rank-i process, across all Relation Algebra (RA) rules
    int *cumulative_tuple_process_map;

    u64 *cumulative_all_to_all_buffer;
    int* cumulative_all_to_all_recv_process_size_array;

    mpi_comm mcomm;


public:

    void push_back(parallel_RA* pj) {RA_list.push_back(pj);}
    void push_relation(relation* G) {relation_manager.push_back(G);}
    void set_refinement_chooser(int rc) {refinement_chooser = rc;}
    void set_refinement_factor(double rf)   {refinement_factor = rf;}
    void set_refinement_interval(int ri)    {refinement_ts = ri;}
    void set_threshold(int thold)   {threshold = thold;}

    void enable_comm_compaction()   {comm_compaction = true;}
    void disable_comm_compaction()  {comm_compaction = false;}
    void set_iteration_count (int ic)   {iteration_count = ic;}
    void allow_print_result ()   {enable_dump_io = true;}
    void push_rule(parallel_RA* pj) {RA_list.push_back(pj);}
    std::vector<relation*> get_relation_manager() {return relation_manager;}
    void set_relation_manager(std::vector<relation*>& rm) {relation_manager = rm;}
    std::vector<parallel_RA*> get_RA_list() {return RA_list;}
    u32 get_bucket_count() {return mcomm.get_local_nprocs();}

    u64 intra_bucket_comm();
    u32 local_compute();
    int all_to_all();
    void local_insert_in_newt();
    u32 local_insert_in_full();
    void check_for_fixed_point(std::vector<u64>& history);
    void print_all_relation();
    void enable_logging();

    void set_comm(mpi_comm& mcomm);
    void execute(int batch_size, std::vector<u64>& history, int task_id);
    void execute_time(double batch_time, std::vector<u64>& history, int task_id);
};

#endif
