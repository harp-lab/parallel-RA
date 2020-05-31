#ifndef LIE_H
#define LIE_H



class LIE
{
private:
    int mode;
    u32 batch_size;
    double batch_time;
    float task_threshold;
    mpi_comm mcomm;
    std::vector<RAM*> tasks;
    std::vector<std::vector<RAM>> taskgraph;

public:
    void push_back(RAM* ra)    {    tasks.push_back(ra);    }

    void set_comm(mpi_comm comm)   { mcomm = comm;  }

    void set_batch_size (u32 bs)    {batch_size = bs;}
    u32 get_batch_size ()    {return batch_size;}

    void set_batch_time (double bt)    {batch_time = bt;}
    double get_batch_time ()    {return batch_time;}

    void set_task_threshold (float th)    {task_threshold = th;}
    float get_task_threshold ()    {return task_threshold;}

    void set_mode (int md)    {mode = md;}
    float get_mode ()    {return mode;}

    void set_taskgraph (std::vector<std::vector<RAM>>& graph)    {taskgraph = graph;}

    bool execute();
};

#endif
