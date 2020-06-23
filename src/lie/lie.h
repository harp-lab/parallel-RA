#ifndef LIE_H
#define LIE_H



class LIE
{
private:
    int mode=1;
    u32 batch_size;
    double batch_time;
    float task_threshold=1.1;
    mpi_comm mcomm;
    std::unordered_set<RAM*> tasks;
    //std::vector<std::vector<RAM>> taskgraph;
    std::unordered_map<RAM*, std::unordered_set<RAM*>> taskgraph1;


public:
    void add_scc(RAM* ra)    {    tasks.insert(ra);    }

    void set_comm(mpi_comm comm)   { mcomm = comm;  }

    void set_batch_size (u32 bs)    {batch_size = bs;}
    u32 get_batch_size ()    {return batch_size;}

    void set_batch_time (double bt)    {batch_time = bt;}
    double get_batch_time ()    {return batch_time;}

    void set_task_threshold (float th)    {task_threshold = th;}
    float get_task_threshold ()    {return task_threshold;}

    void set_mode (int md)    {mode = md;}
    float get_mode ()    {return mode;}

    //void set_taskgraph (std::vector<std::vector<RAM>>& graph)    {taskgraph = graph;}

    void add_scc_dependance (RAM* src_task, RAM* destination_task);

    bool execute();
};

#endif
