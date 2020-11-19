/*
 * scc (tasks)
 * Copyright (c) Sidharth Kumar, et al, see License.md
 */



#ifndef RAM_H
#define RAM_H



class RAM
{

private:

    int ram_id;
    bool init_status=false;

    int iteration_count = -1;                               /// Number of iterations in a fixed point

    bool logging = false;                                   /// If logging is enabled or not

    double refinement_factor = 4;                           /// For spatial balancing, a bucket is broken into 4 sub-buckets


    u32 ram_relation_count;
    relation *ram_relations[64];
    bool ram_relation_status[64];


    std::vector<parallel_RA*> RA_list;                      /// All relations of this SCC

    u64 *intra_bucket_buf_output_size;                      /// results of intra-bucket comm
    u64 **intra_bucket_buf_output;

    all_to_allv_buffer compute_buffer;                       /// result of compute

    all_to_all_buffer all_to_all_compute_buffer;                       /// result of compute

    u64 *cumulative_all_to_allv_buffer;                      /// result of all to all comm
    int* cumulative_all_to_allv_recv_process_size_array;

    mpi_comm mcomm;                                         /// comm related

    u32 loop_count_tracker;

public:

    ~RAM();
    RAM (bool ic, int ram_id);



    /// Set local task-level communicator
    void set_comm(mpi_comm& mcomm);


    /// For debugging purpose
    int get_id() {return ram_id;}


    /// add relations pertaining to this SCC
    /// void add_relation(relation* G) {ram_relations.insert(G);}



    //std::map<u32, std::map<relation*, bool>> get_RAM_relations()   {return ram_relations;}
    relation** get_RAM_relations()   {return ram_relations;}
    bool* get_RAM_relations_status()   {return ram_relation_status;}


    /// add relations pertaining to this SCC
    void add_relation(relation*& G, bool i_status);


    /// add rule to the SCC
    void add_rule(parallel_RA* pj) {RA_list.push_back(pj);}


    /// Load balancing related
    void set_refinement_factor(double rf)   {refinement_factor = rf;}


    /// Iteration count set through the constructor
    void set_iteration_count (int ic)   {iteration_count = ic;}
    int get_iteration_count ()   {return iteration_count;}


    /// for debugging
    void enable_logging()   {logging = true;}
    void print_all_relation();


    void io_all_relation(int status);


    /// the buckets over which the SCC is spread across
    u32 get_bucket_count() {return mcomm.get_local_nprocs();}


    u32 get_ram_relation_count() {return ram_relation_count;}


    /// Spatial balancing of all relations
    void load_balance();


    /// Intra bucket comm for sub-buckets
    u64 intra_bucket_comm_execute();


    /// Buffer to hold new tuples
    u32 allocate_compute_buffers();


    void allocate_all_to_all_compute_buffers();


    u32 local_compute_with_all_to_all_threshold(int* offset);

    /// Join/compy/acopy
    u32 local_compute(int* offset);


    //// Comm compaction
    void all_to_all();
    void all_to_all_with_threshold();


    /// Free intermediate buffers
    void free_compute_buffers();
    void free_all_to_all_compute_buffers();


    /// Update the head relation with new tuples
    void local_insert_in_newt(std::map<u64, u64>& intern_map);
    void local_insert_in_newt_with_threshold(std::map<u64, u64>& intern_map);


    /// insert delta in full, copy newt pointer to delta
    void local_insert_in_full();


    void insert_delta_in_full();

    /// has fixed point reached
    void check_for_fixed_point(std::vector<u32>& history);


    void execute_in_batches_with_threshold(int batch_size, std::vector<u32>& history, std::map<u64, u64>& intern_map, double *running_time);


    /// Start running this SCC (task) for "batck_size" iterations
    void execute_in_batches(int batch_size, std::vector<u32>& history, std::map<u64, u64>& intern_map, double *running_time, double *running_intra_bucket_comm, double *running_buffer_allocate, double *running_local_compute, double *running_all_to_all, double *running_buffer_free, double *running_insert_newt, double *running_insert_in_full, double *running_fp);



    void execute_in_batches_with_all_to_all_threshold(int batch_size, std::vector<u32>& history, std::map<u64, u64>& intern_map, double *running_time, double *running_intra_bucket_comm, double *running_buffer_allocate, double *running_local_compute, double *running_all_to_all, double *running_buffer_free, double *running_insert_newt, double *running_insert_in_full);


    /// Start running this SCC (task) for "batck_time" seconds
    void execute_by_wall_clock(double batch_time, std::vector<u32>& history, std::map<u64, u64>& intern_map, double *running_time, double *running_intra_bucket_comm, double *running_buffer_allocate, double *running_local_compute, double *running_all_to_all, double *running_buffer_free, double *running_insert_newt, double *running_insert_in_full);
};

#endif
