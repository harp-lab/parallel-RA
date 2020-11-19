/*
 * Logical Inferencing Engine (LIE)
 * Copyright (c) Sidharth Kumar, et al, see License.md
 */



#ifndef LIE_H
#define LIE_H



class LIE
{
private:

    u32 batch_size=10;                                                 /// Number of iterations between two fixed-point checks
    double batch_time;                                              /// Wallclock time between two fixed-point checks
    mpi_comm mcomm;                                                 /// MPI class

    u32 lie_relation_count;
    relation *lie_relations[64];

    //u32 lie_relations_key;
    //std::map<u32, relation*> lie_relations;                    /// List of all relations

    u32 lie_sccs_count;
    RAM *lie_sccs[64];

    //u32 lie_sccs_key;
    //std::map<u32, RAM*> lie_sccs;                                 /// List of all tasks

    std::map<RAM*, std::set<RAM*>> taskgraph;  /// List of all edges in a task (this is an adjacency list)

    std::map<u64, u64> intern_map;                        /// Intern table

    bool enable_io;

    bool restart_flag;

    std::string restart_dir_name;                          /// the directory of restart

    bool share_io;                                        /// whether using MPI collective IO to write file

    bool offset_io;										 /// whether read checkpoint dump with offset

    bool separate_io;                                    /// whether write checkpoint dump separately for each process

    int loop_counter;

    std::vector<int> executed_scc_id;

    std::string output_dir;

    int cp_iteration;

public:

    ~LIE();

    LIE()
    {
    	cp_iteration = 0;
    	executed_scc_id = {};
    	loop_counter = 0;
    	separate_io = false;
    	offset_io = false;
    	share_io = false;
    	restart_flag = false;
        enable_io = false;
        lie_relation_count = 0;
        lie_sccs_count = 0;
        taskgraph = {{},{}};
        intern_map = {{},{}};
    }

    void set_cp_iteration(int iteration)       {cp_iteration = iteration;}

    void set_output_dir(std::string output)   {output_dir = output;}

    void set_executed_scc_id(std::vector<int> id)     {executed_scc_id = id;}

    std::vector<int> get_executed_scc_id()       {return executed_scc_id;}

    void set_loop_counter(int count)      {loop_counter = count;}

    int get_loop_counter()         {return loop_counter;}

    void enable_separate_io()      {separate_io = true;}

    void enable_offset_io()    {offset_io = true;}

    void enable_share_io()    {share_io = true;}

    void set_restart_dir_name(std::string path)  {restart_dir_name = path;}

    void set_restart_flag(bool flag)   {restart_flag = flag;}

    void enable_IO()    {enable_io = true;}


    void print_all_relation_size();

    /// Sets the communicator object
    void set_comm(mpi_comm comm)   { mcomm = comm;  }


    /// Batch size
    void set_batch_size (u32 bs)    {batch_size = bs;}
    u32 get_batch_size ()    {return batch_size;}


    /// Batch time
    void set_batch_time (double bt)    {batch_time = bt;}
    double get_batch_time ()    {return batch_time;}


    /// Adds a new relation to the LIE
    void add_relation(relation* rel);


    /// Adds a new SCC to the LIE
    void add_scc(RAM* ra);


    /// Returns a list of infiished tasks
    RAM* one_runnable_tasks();


    /// Removes tasks from the task edge list and tas list
    void update_task_graph(RAM* removable_tasks);


    /// Populates the adjacency list of tasks
    void add_scc_dependance (RAM* src_task, RAM* destination_task);


    /// Runs all tasks within the LIE, following the dependence as set by taskgraph
    bool execute();

    void write_checkpoint_dump(int loop_counter, std::vector<int> executed_scc_id);
};

#endif
