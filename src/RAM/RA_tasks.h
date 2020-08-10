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

    std::string name;                                       /// Name of this task (for debugging only)
    int iteration_count = -1;                               /// Number of iterations in a fixed point

    bool logging = false;                                   /// If logging is enabled or not

    double refinement_factor = 4;                           /// For spatial balancing, a bucket is broken into 4 sub-buckets


    std::unordered_map<relation*, bool> ram_relations;            /// All relations of this SCC
    std::vector<parallel_RA*> RA_list;                      /// All relations of this SCC

    u64 *intra_bucket_buf_output_size;                      /// results of intra-bucket comm
    u64 **intra_bucket_buf_output;

    all_to_all_buffer compute_buffer;                       /// result of compute

    u64 *cumulative_all_to_all_buffer;                      /// result of all to all comm
    int* cumulative_all_to_all_recv_process_size_array;

    mpi_comm mcomm;                                         /// comm related

    u32 loop_count_tracker;

public:

    ~RAM();
    RAM()   {}
    RAM (bool ic, int ram_id);
    RAM (bool ic);
    RAM (bool ic, std::string name);


    /// Set local task-level communicator
    void set_comm(mpi_comm& mcomm);


    /// For debugging purpose
    std::string get_name() {return name;}


    /// For debugging purpose
    int get_id() {return ram_id;}


    /// add relations pertaining to this SCC
    /// void add_relation(relation* G) {ram_relations.insert(G);}



    std::unordered_map<relation*, bool> get_RAM_relations()   {return ram_relations;}


    /// add relations pertaining to this SCC
    void add_relation(relation* G, bool i_status)
    {
        //init_status = i_status;
        //ram_relations.insert(G);
        //ram_relations[G] = i_status;
        auto it = ram_relations.find(G);
        if( it != ram_relations.end() ) {
            it->second = i_status;
        }
        else {
            ram_relations.insert(std::make_pair(G, i_status));
        }

    }


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


    /// the buckets over which the SCC is spread across
    u32 get_bucket_count() {return mcomm.get_local_nprocs();}


    /// Spatial balancing of all relations
    void load_balance();


    /// Intra bucket comm for sub-buckets
    u64 intra_bucket_comm_execute();


    /// Buffer to hold new tuples
    u32 allocate_compute_buffers();


    /// Join/compy/acopy
    u32 local_compute();


    //// Comm compaction
    void all_to_all();


    /// Free intermediate buffers
    void free_compute_buffers();


    /// Update the head relation with new tuples
    void local_insert_in_newt(std::unordered_map<u64, u64>& intern_map);


    /// insert delta in full, copy newt pointer to delta
    void local_insert_in_full();


    void insert_delta_in_full();

    /// has fixed point reached
    void check_for_fixed_point(std::vector<u32>& history);


    /// Start running this SCC (task) for "batck_size" iterations
    void execute_in_batches(int batch_size, std::vector<u32>& history, std::unordered_map<u64, u64>& intern_map);


    /// Start running this SCC (task) for "batck_time" seconds
    void execute_by_wall_clock(double batch_time, std::vector<u32>& history, std::unordered_map<u64, u64>& intern_map);
};

#endif
