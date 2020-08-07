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

    std::unordered_set<relation*> lie_relations;                    /// List of all relations
    std::unordered_set<RAM*> tasks;                                 /// List of all tasks
    std::unordered_map<RAM*, std::unordered_set<RAM*>> taskgraph1;  /// List of all edges in a task (this is an adjacency list)

    std::unordered_map<u64, u64> intern_map;                        /// Intern table


public:

    void print_all_relation();

    /// Sets the communicator object
    void set_comm(mpi_comm comm)   { mcomm = comm;  }


    /// Batch size
    void set_batch_size (u32 bs)    {batch_size = bs;}
    u32 get_batch_size ()    {return batch_size;}


    /// Batch time
    void set_batch_time (double bt)    {batch_time = bt;}
    double get_batch_time ()    {return batch_time;}


    /// Adds a new relation to the LIE
    void add_relation(relation* rel)    {    lie_relations.insert(rel);    }


    /// Adds a new SCC to the LIE
    void add_scc(RAM* ra)    {    tasks.insert(ra);    }


    /// Returns a list of infiished tasks
    std::unordered_set<RAM*> list_of_runnable_tasks(std::unordered_set<RAM*> tasks, std::unordered_map<RAM*, std::unordered_set<RAM*>> taskgraph1);


    /// Removes tasks from the task edge list and tas list
    void update_task_graph(std::unordered_set<RAM*> removable_tasks);


    /// Populates the adjacency list of tasks
    void add_scc_dependance (RAM* src_task, RAM* destination_task);


    /// Runs all tasks within the LIE, following the dependence as set by taskgraph1
    bool execute();
};

#endif
