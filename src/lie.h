#ifndef LIE_H
#define LIE_H

#include <mpi.h>
#include <vector>
#include <unordered_set>
#include <queue>
#include "RA_tasks.h"


class LIE
{
private:
    int mode;
    u32 batch_size;
    float task_threshold;
    mpi_comm mcomm;
    std::vector<RAM*> tasks;
    std::vector<std::vector<RAM>> taskgraph;

public:
    void push_back(RAM* ra)    {    tasks.push_back(ra);    }

    void set_comm(mpi_comm comm)   { mcomm = comm;  }

    void set_batch_size (u32 bs)    {batch_size = bs;}
    u32 get_batch_size ()    {return batch_size;}
    void set_task_threshold (float th)    {task_threshold = th;}
    float get_task_threshold ()    {return task_threshold;}

    void set_mode (int md)    {mode = md;}
    float get_mode ()    {return mode;}

    void set_taskgraph (std::vector<std::vector<RAM>>& graph)    {taskgraph = graph;}

    bool execute ()
    {
        double init_start = MPI_Wtime();
        int color = -1;
        u32 stratified_task_layers = taskgraph.size();

        int number_of_parallel_tasks = taskgraph[0].size();
        u64 *tuple_count_per_task = new u64[number_of_parallel_tasks];
        int *ranks_per_task = new int[number_of_parallel_tasks];
        int *cumulative_ranks_per_task = new int[number_of_parallel_tasks + 1];
        memset(cumulative_ranks_per_task, 0, (number_of_parallel_tasks + 1) * sizeof(int));

        MPI_Comm newcomm;
        u64 total_tuple_count = 0;

        // Initialization strata
        for (int j=0; j < number_of_parallel_tasks; j++)
        {
            RAM task = taskgraph[0][j];
            tuple_count_per_task[j] = (get_total_number_of_tuples(task));
            total_tuple_count = total_tuple_count + tuple_count_per_task[j];
        }
#if 1
        u32 temp_rank_count1 = 0;
        u32 temp_rank_count2 = 0;
        for (int j=0; j < number_of_parallel_tasks; j++)
            temp_rank_count1 = temp_rank_count1 + ceil(mcomm.get_nprocs() * (float)((float)tuple_count_per_task[j]/total_tuple_count));
        u32 deficiet = temp_rank_count1 - mcomm.get_nprocs();

        cumulative_ranks_per_task[0] = 0;
        for (int j=0; j < number_of_parallel_tasks; j++)
        {
            if (tuple_count_per_task[j] == 0)
            {
                ranks_per_task[j] = 0;
                cumulative_ranks_per_task[j+1] = (cumulative_ranks_per_task[j]);
            }
            else if (deficiet != 0)
            {
                int temp = ceil(mcomm.get_nprocs() * (float)((float)tuple_count_per_task[j]/total_tuple_count));
                if (temp == 1)
                {
                    ranks_per_task[j] = (ceil(mcomm.get_nprocs() * (float)((float)tuple_count_per_task[j]/total_tuple_count)));
                    cumulative_ranks_per_task[j+1] = (cumulative_ranks_per_task[j] + ranks_per_task[j]);
                    temp_rank_count2 = temp_rank_count2 + ranks_per_task[j];
                }
                else
                {
                    ranks_per_task[j] = (floor(mcomm.get_nprocs() * (float)((float)tuple_count_per_task[j]/total_tuple_count)));
                    cumulative_ranks_per_task[j+1] = (cumulative_ranks_per_task[j] + ranks_per_task[j]);
                    temp_rank_count2 = temp_rank_count2 + ranks_per_task[j];
                    deficiet--;
                }
            }
            else
            {
                ranks_per_task[j] = (ceil(mcomm.get_nprocs() * (float)((float)tuple_count_per_task[j]/total_tuple_count)));
                cumulative_ranks_per_task[j+1] = (cumulative_ranks_per_task[j] + ranks_per_task[j]);
                temp_rank_count2 = temp_rank_count2 + ranks_per_task[j];
            }
        }

#else
        cumulative_ranks_per_task[0] = 0;
        for (int j=0; j < number_of_parallel_tasks; j++)
        {
            ranks_per_task[j] = (ceil(mcomm.get_nprocs() * (float)((float)tuple_count_per_task[j]/total_tuple_count)));
            cumulative_ranks_per_task[j+1] = (cumulative_ranks_per_task[j] + ranks_per_task[j]);

            if (cumulative_ranks_per_task[j+1] > mcomm.get_nprocs())
            {
                ranks_per_task[j] = mcomm.get_nprocs() - cumulative_ranks_per_task[j];
                cumulative_ranks_per_task[j+1] = cumulative_ranks_per_task[j] + ranks_per_task[j];
                assert(j == number_of_parallel_tasks-1);
            }
        }
#endif


        for (int j=0; j < number_of_parallel_tasks; j++)
        {
            if (mcomm.get_rank() < cumulative_ranks_per_task[j+1])
            {
                color = j;
                break;
            }
        }


        MPI_Comm_split(mcomm.get_comm(), color, mcomm.get_rank(), &newcomm);
        mcomm.set_local_comm(&newcomm);

        for (u32 i=0; i < stratified_task_layers; i++)
            for (int j=0; j < number_of_parallel_tasks; j++)
                taskgraph[i][j].set_comm(mcomm);


        initialize_relations(taskgraph[0][color]);

#if 1
        int loop_count = 0;

        double exec_time = 0;
        double rebalance_check_time = 0;
        double rebalance_time = 0;
        double init_end = MPI_Wtime();

        if (mcomm.get_rank() == 0)
            std::cout << "Initialization Time: " << (init_end - init_start) << std::endl;

        for (u32 i=0; i < stratified_task_layers; i++)
        {
            //int finished_task_count = 0;
            mpi_comm new_mcomm(mcomm);

            std::vector<u64> history;
            int finished_task_count=0;

loop_again:

            double loop_exec_start = MPI_Wtime();
            number_of_parallel_tasks = taskgraph[i].size();
            taskgraph[i][color].execute(batch_size, history);
            double loop_exec_end = MPI_Wtime();
            exec_time = exec_time + (loop_exec_end - loop_exec_start);

            if (mcomm.get_rank() == 0)
                std::cout << "OUTER E [" << color << "] Loop " << loop_count << " Exec time: " << (loop_exec_end - loop_exec_start) << " Total exec time " << exec_time << std::endl;

#if 0
            if (mcomm.get_rank() == 0)
            {
                std::cout << "History ";
                for (u32 x=0; x < history.size(); x++)
                    std::cout << history[x] << " ";
                std::cout << std::endl;
            }
#endif


            if (i == 0) continue;

#if 1
            double rebalance_check_start = MPI_Wtime();
            int new_color = -1;
            int* prev_ranks_per_task = new int[number_of_parallel_tasks];
            memcpy(prev_ranks_per_task, ranks_per_task, number_of_parallel_tasks * sizeof(int));
            bool check_for_rebalance = rebalance_comm(history, number_of_parallel_tasks, tuple_count_per_task, ranks_per_task, cumulative_ranks_per_task, new_mcomm, &color, &new_color, &finished_task_count, mode);
            double rebalance_check_end = MPI_Wtime();
            rebalance_check_time = rebalance_check_time + (rebalance_check_end - rebalance_check_start);

            if (mcomm.get_rank() == 0)
                std::cout << "OUTER RC [" << color << "] Loop " << loop_count << " Rebalance check time: " << (rebalance_check_end - rebalance_check_start) << " Total exec time " << rebalance_check_time << std::endl;


            double rebalance_start = MPI_Wtime();
            if (check_for_rebalance == true)
            {
                rebalance_data(taskgraph[i][color], taskgraph[i][new_color], cumulative_ranks_per_task[color], new_mcomm, tuple_count_per_task[color], prev_ranks_per_task[color], ranks_per_task[color]);
                taskgraph[i][new_color].set_comm(new_mcomm);
                color = new_color;
            }
            delete[] prev_ranks_per_task;
            double rebalance_end = MPI_Wtime();
            rebalance_time = rebalance_time + (rebalance_end - rebalance_start);
#endif

            if (mcomm.get_rank() == 0)
                std::cout << "OUTER R [" << color << "] Loop " << loop_count
                          << " Finished task count: " << finished_task_count
                          << " Rebalance time: " << (rebalance_end - rebalance_start)
                          << " Total exec time " << rebalance_time
                          << std::endl << std::endl;

            if (finished_task_count != number_of_parallel_tasks)
            {
                loop_count++;
                //if (mcomm.get_rank() == 0)
                //    std::cout << "Loop count " << loop_count << " Finished task count " << finished_task_count << " Number of parallel tasks " << number_of_parallel_tasks << std::endl;
                goto loop_again;
            }
        }
#endif
        double finalize_start = MPI_Wtime();
        finalize_relation(taskgraph[0][color]);

        delete[] tuple_count_per_task;
        delete[] ranks_per_task;
        delete[] cumulative_ranks_per_task;


        double finalize_end = MPI_Wtime();

        if (mcomm.get_rank() == 0)
            std::cout << "Finalize Time: " << (finalize_end - finalize_start) << std::endl << std::endl;

        double total_time = finalize_end - init_start;
        double max_time = 0;
        MPI_Allreduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, mcomm.get_comm());

        if (total_time == max_time)
        {
            std::cout << "[M] [" << mcomm.get_rank() << "] Total Time: " << (finalize_end - finalize_start) + (rebalance_time) + rebalance_check_time + exec_time + (init_end - init_start)
                      << " Init time: " << (init_end - init_start)
                      << " Exec time: " << exec_time
                      << " Rebalance check time: " << rebalance_check_time
                      << " Rebalance time: " << rebalance_time
                      << std::endl;
        }

        MPI_Barrier(MPI_COMM_WORLD);
        //if (mcomm.get_rank() == 0)
        {
            std::cout << "[" << mcomm.get_rank() << "] Total Time: " << (finalize_end - finalize_start) + (rebalance_time) + rebalance_check_time + exec_time + (init_end - init_start)
                      << " Init time: " << (init_end - init_start)
                      << " Exec time: " << exec_time
                      << " Rebalance check time: " << rebalance_check_time
                      << " Rebalance time: " << rebalance_time
                      << std::endl;
        }

        return true;
    }

    void rebalance_data(RAM& send_task, RAM& recv_task, int rank_offset, mpi_comm new_mcomm, int tuples_per_task, int ib, int ob)
    {
        std::vector<relation*> send_rm = send_task.get_relation_manager();
        std::vector<relation*> recv_rm = recv_task.get_relation_manager();
        for (u32 r = 0; r < send_rm.size(); r++)
        {
            relation* send_rel =  send_rm[r];
            relation* recv_rel =  recv_rm[r];
            send_rel->copy_relation(recv_rel, new_mcomm, rank_offset, tuples_per_task, ib, ob);
        }
    }


    bool rebalance_comm(std::vector<u64>& history, int number_of_parallel_tasks, u64*& current_tuple_count_per_task, int*& current_ranks_per_task, int*& current_cumulative_ranks_per_task, mpi_comm& new_mcomm, int *current_color, int *next_color, int* finished_task_count, int mode)
    {

        // Check if balancing is required

        u64 *next_tuple_count_per_task_local = new u64[number_of_parallel_tasks];
        u64 *next_tuple_count_per_task = new u64[number_of_parallel_tasks];
        memset(next_tuple_count_per_task_local, 0, number_of_parallel_tasks * sizeof(u64));
        memset(next_tuple_count_per_task, 0, number_of_parallel_tasks * sizeof(u64));

        // D F D F D F D F D F
        // 0 1 2 3 4 5 6 7 8 9
        u64 full_in_scc = 0;
        u64 delta_in_scc = 0;

        if (mode == -1)
        {
            full_in_scc = history[history.size()-1];
            delta_in_scc = history[history.size()-2];
            if (delta_in_scc == 0)  full_in_scc = 0;

            next_tuple_count_per_task_local[*current_color] = full_in_scc;
            MPI_Allreduce(next_tuple_count_per_task_local, next_tuple_count_per_task, number_of_parallel_tasks, MPI_UNSIGNED_LONG_LONG, MPI_BOR, mcomm.get_comm());
            memset(current_tuple_count_per_task, 0, number_of_parallel_tasks * sizeof(u64));

            for (int j=0; j < number_of_parallel_tasks; j++)
                current_tuple_count_per_task[j] = next_tuple_count_per_task[j];

            delete[] next_tuple_count_per_task_local;
            delete[] next_tuple_count_per_task;
            return false;
        }


        if (history.size() > 6 && mode == 1)
        {
            full_in_scc = history[history.size()-1];
            delta_in_scc = history[history.size()-2];
            if (delta_in_scc == 0)  full_in_scc = 0;

            next_tuple_count_per_task_local[*current_color] = full_in_scc;
            MPI_Allreduce(next_tuple_count_per_task_local, next_tuple_count_per_task, number_of_parallel_tasks, MPI_UNSIGNED_LONG_LONG, MPI_BOR, mcomm.get_comm());

            // x = 2c - p
            for (int j=0; j < number_of_parallel_tasks; j++)
            {
                if (current_tuple_count_per_task[j] >=  2 * next_tuple_count_per_task[j])
                    next_tuple_count_per_task[j] = 0;
                else
                    next_tuple_count_per_task[j] = 2 * next_tuple_count_per_task[j] - current_tuple_count_per_task[j];
            }
        }
        else
        {
            full_in_scc = history[history.size()-1];
            delta_in_scc = history[history.size()-2];
            if (delta_in_scc == 0)  full_in_scc = 0;

            next_tuple_count_per_task_local[*current_color] = full_in_scc;
            MPI_Allreduce(next_tuple_count_per_task_local, next_tuple_count_per_task, number_of_parallel_tasks, MPI_UNSIGNED_LONG_LONG, MPI_BOR, mcomm.get_comm());
        }



        delete[] next_tuple_count_per_task_local;


        u64 total_tuple_count = 0;

        (*finished_task_count) = 0;
        for (int j=0; j < number_of_parallel_tasks; j++)
        {
            total_tuple_count = total_tuple_count + next_tuple_count_per_task[j];
            if (next_tuple_count_per_task[j] == 0)
                (*finished_task_count)++;
        }

        int* next_cumulative_ranks_per_task = new int[number_of_parallel_tasks + 1];
        memset(next_cumulative_ranks_per_task, 0, (number_of_parallel_tasks + 1) * sizeof(int));

        int* next_ranks_per_task = new int[number_of_parallel_tasks];



#if 1
        u32 temp_rank_count1 = 0;
        u32 temp_rank_count2 = 0;
        for (int j=0; j < number_of_parallel_tasks; j++)
            temp_rank_count1 = temp_rank_count1 + ceil(mcomm.get_nprocs() * (float)((float)next_tuple_count_per_task[j]/total_tuple_count));
        u32 deficiet = temp_rank_count1 - mcomm.get_nprocs();

        next_cumulative_ranks_per_task[0] = 0;
        for (int j=0; j < number_of_parallel_tasks; j++)
        {
            if (next_tuple_count_per_task[j] == 0)
            {
                next_ranks_per_task[j] = 0;
                next_cumulative_ranks_per_task[j+1] = (next_cumulative_ranks_per_task[j]);
            }
            else if (deficiet != 0)
            {
                int temp = ceil(mcomm.get_nprocs() * (float)((float)next_tuple_count_per_task[j]/total_tuple_count));
                if (temp == 1)
                {
                    next_ranks_per_task[j] = (ceil(mcomm.get_nprocs() * (float)((float)next_tuple_count_per_task[j]/total_tuple_count)));
                    next_cumulative_ranks_per_task[j+1] = (next_cumulative_ranks_per_task[j] + next_ranks_per_task[j]);
                    temp_rank_count2 = temp_rank_count2 + next_ranks_per_task[j];
                }
                else
                {
                    next_ranks_per_task[j] = (floor(mcomm.get_nprocs() * (float)((float)next_tuple_count_per_task[j]/total_tuple_count)));
                    next_cumulative_ranks_per_task[j+1] = (next_cumulative_ranks_per_task[j] + next_ranks_per_task[j]);
                    temp_rank_count2 = temp_rank_count2 + next_ranks_per_task[j];
                    deficiet--;
                }
            }
            else
            {
                next_ranks_per_task[j] = (ceil(mcomm.get_nprocs() * (float)((float)next_tuple_count_per_task[j]/total_tuple_count)));
                next_cumulative_ranks_per_task[j+1] = (next_cumulative_ranks_per_task[j] + next_ranks_per_task[j]);
                temp_rank_count2 = temp_rank_count2 + next_ranks_per_task[j];
            }
        }
#else
        next_cumulative_ranks_per_task[0] = 0;
        for (int j=0; j < number_of_parallel_tasks; j++)
        {
            next_ranks_per_task[j] = (ceil(mcomm.get_nprocs() * (float)((float)next_tuple_count_per_task[j]/total_tuple_count)));
            next_cumulative_ranks_per_task[j+1] = (next_cumulative_ranks_per_task[j] + next_ranks_per_task[j]);

            if (next_cumulative_ranks_per_task[j+1] > mcomm.get_nprocs())
            {
                next_ranks_per_task[j] = mcomm.get_nprocs() - next_cumulative_ranks_per_task[j];
                next_cumulative_ranks_per_task[j+1] = next_cumulative_ranks_per_task[j] + next_ranks_per_task[j];

                if (j != number_of_parallel_tasks-1) std::cout << "J is " << j << std::endl;
                assert(j == number_of_parallel_tasks-1);
            }
        }
#endif

        int count=0;
        for (int j=0; j < number_of_parallel_tasks; j++)
            if (current_ranks_per_task[j] == next_ranks_per_task[j])
                count++;

        float max = -10.0f;
        int num, den;
        for (int i=0; i < number_of_parallel_tasks; i++)
        {
            if (next_ranks_per_task[i] == 0 || current_ranks_per_task[i] == 0)
                continue;

            if (max < (float)((float)current_ranks_per_task[i]/(float)next_ranks_per_task[i]))
            {
                max = (float)((float)current_ranks_per_task[i]/(float)next_ranks_per_task[i]);
                num = current_ranks_per_task[i];
                den = next_ranks_per_task[i];
            }

            if (max < (float)((float)next_ranks_per_task[i]/(float)current_ranks_per_task[i]))
            {
                max = (float)((float)next_ranks_per_task[i]/(float)current_ranks_per_task[i]);
                den = current_ranks_per_task[i];
                num = next_ranks_per_task[i];
            }
        }





        if (total_tuple_count == 0)
        {
            if (mcomm.get_rank() == 0)
                std::cout << "-------------------All Tasks finished-------------------" << std::endl;

            delete[] next_ranks_per_task;
            delete[] next_cumulative_ranks_per_task;
            return false;
        }

        if (count == number_of_parallel_tasks)
        {
            if (mcomm.get_rank() == 0)
            {
                std::cout << "-------------------No change in rank distribution-------------------" << std::endl;
                for (int j=0; j < number_of_parallel_tasks; j++)
                    std::cout << "[" << current_ranks_per_task[j] << " " << next_ranks_per_task[j] << "] ";
                std::cout << std::endl;
            }

            delete[] next_ranks_per_task;
            delete[] next_cumulative_ranks_per_task;
            return false;
        }

        if (max < task_threshold)
        {
            if (mcomm.get_rank() == 0)
            {
                std::cout << "Does not cross threshold of 1.1 current threshold is = " << max << " " << num << " " << den << std::endl;
                for (int j=0; j < number_of_parallel_tasks; j++)
                    std::cout << "[" << current_ranks_per_task[j] << " " << next_ranks_per_task[j] << "] ";
                std::cout << std::endl;
            }
            delete[] next_cumulative_ranks_per_task;
            return false;
        }
        else
        {
            if (mcomm.get_rank() == 0){
                std::cout << "-----------------Crosses the threshold of 1.1 current threshold is = " << max << " " << num << " " << den << "-----------------" << std::endl;
                for (int j=0; j < number_of_parallel_tasks; j++)
                    std::cout << "[" << current_ranks_per_task[j] << " " << next_ranks_per_task[j] << "] ";
                std::cout << std::endl;
            }
        }


        *next_color = -1;
        for (int j=0; j < number_of_parallel_tasks; j++)
        {
            if (mcomm.get_rank() < next_cumulative_ranks_per_task[j+1])
            {
                *next_color = j;
                break;
            }
        }


        MPI_Comm newcomm;
        MPI_Comm_split(mcomm.get_comm(), *next_color, mcomm.get_rank(), &newcomm);
        new_mcomm.set_local_comm(&newcomm);

        for (int j=0; j < number_of_parallel_tasks; j++)
        {
            current_ranks_per_task[j] = next_ranks_per_task[j];
            current_tuple_count_per_task[j] = next_tuple_count_per_task[j];
            current_cumulative_ranks_per_task[j] = next_cumulative_ranks_per_task[j];
        }
        current_cumulative_ranks_per_task[number_of_parallel_tasks] = next_cumulative_ranks_per_task[number_of_parallel_tasks];

        delete[] next_ranks_per_task;
        delete[] next_cumulative_ranks_per_task;
        delete[] next_tuple_count_per_task;

        return true;
    }


    void initialize_relations(RAM& task)
    {
        std::vector<relation*> rm = task.get_relation_manager();
        for (relation* rel : rm)
            rel->initialize_relation(mcomm);
    }


    void finalize_relation(RAM& task)
    {
        std::vector<relation*> rm = task.get_relation_manager();
        for (relation* rel : rm)
            rel->finalize_relation();
    }


    u64 get_total_number_of_tuples(RAM& task)
    {
        u64 tuple_count = 0;
        std::vector<relation*> rm = task.get_relation_manager();
        for (relation* rel : rm)
        {
            tuple_count = tuple_count + rel->get_delta_element_count();
            tuple_count = tuple_count + rel->get_full_element_count();
        }
        return tuple_count;
    }


    u64 get_total_number_of_tuples_per_task(RAM& task)
    {
        u64 tuple_count = 0;
        u64 global_tuple_count = 0;
        std::vector<relation*> rm = task.get_relation_manager();
        for (relation* rel : rm)
        {
            tuple_count = tuple_count + rel->get_delta_element_count();
            tuple_count = tuple_count + rel->get_full_element_count();
        }

        MPI_Allreduce(&tuple_count, &global_tuple_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mcomm.get_local_comm());
        return global_tuple_count;

    }
};

#endif
