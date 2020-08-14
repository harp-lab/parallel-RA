/*
 * Logical Inferencing Engine (LIE)
 * Copyright (c) Sidharth Kumar, et al, see License.md
 */



#ifndef LIE_H
#define LIE_H



class LIE
{
private:

    u32 batch_size;                                                 /// Number of iterations between two fixed-point checks
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


public:

    ~LIE();

    LIE()
    {
        lie_relation_count = 0;
        //lie_relations_key = 0;
        //lie_relations = {{},{}};

        lie_sccs_count = 0;
        //lie_sccs_key = 0;
        //lie_sccs = {{},{}};

        taskgraph = {{},{}};
        intern_map = {{},{}};
    }

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
};

#endif
