/*
 * Logical Inferencing Engine (LIE)
 * Copyright (c) Sidharth Kumar, et al, see License.md
 */



#include "../parallel_RA_inc.h"


/// This function currentl only returns one runnable task
/// size of returned list is always going to be 1
std::unordered_set<RAM*> LIE::list_of_runnable_tasks(std::unordered_set<RAM*> tasks, std::unordered_map<RAM*, std::unordered_set<RAM*>> taskgraph1)
{
    std::unordered_set<RAM*> runnable_tasks;
    for ( auto itg = tasks.begin(); itg != tasks.end(); ++itg )
    {
        RAM* task = (*itg);
        bool break_loop = false;

        for (auto it = taskgraph1.begin(); it != taskgraph1.end(); it++)
        {
            std::unordered_set<RAM*> it2 = it->second;
            for (auto dit2 = it2.begin(); dit2 != it2.end(); dit2++)
            {
                if (*itg == *dit2)
                {
                    break_loop=true;
                    break;
                }
            }
            if (break_loop==true)
                break;
        }
        if (break_loop==false)
            runnable_tasks.insert(task);

        /// returns as soon as one runnable task is found
        if (runnable_tasks.size() == 1)
            break;
    }
    return runnable_tasks;
}



void LIE::update_task_graph(std::unordered_set<RAM*> executable_task)
{
    for ( auto itg = executable_task.begin(); itg != executable_task.end(); ++itg )
    {
        tasks.erase(*itg);
        taskgraph1.erase(*itg);
    }
}



void LIE::add_scc_dependance (RAM* src_task, RAM* destination_task)
{
    auto it = taskgraph1.find(src_task);
    if( it != taskgraph1.end() )
    {
        auto it2 = (it->second).find(destination_task);
        if( it2 == (it->second).end() )
        {
            (it->second).insert(destination_task);
            taskgraph1[src_task] = it->second;
        }
    }
    else
    {
        std::unordered_set<RAM*> k;
        k.insert(destination_task);
        taskgraph1.insert(std::make_pair(src_task, k));
    }

    std::unordered_set<RAM*> temp = taskgraph1[src_task];
    temp.insert(destination_task);

    return;
}


void LIE::print_all_relation()
{
    u64 total_facts=0;
    for (std::unordered_set<relation*>::iterator it = lie_relations.begin() ; it != lie_relations.end(); ++it)
    {
        std::cout << (*it)->get_debug_id() << ": {" << (*it)->get_arity() << "}. (" << (*it)->get_full_element_count() << " total facts)" << std::endl;
        total_facts = total_facts + (*it)->get_full_element_count();
    }
    std::cout << "Total facts across all relations " << total_facts << std::endl;
}


bool LIE::execute ()
{

    /// Main : Execute : init : start

    mcomm.set_local_comm(mcomm.get_comm());


    /// Initialize all relations
    for (std::unordered_set<relation*>::iterator it = lie_relations.begin() ; it != lie_relations.end(); ++it)
        (*it)->initialize_relation(mcomm);


    /// Executable task
    std::unordered_set<RAM*> executable_task = list_of_runnable_tasks(tasks, taskgraph1);
    assert(executable_task.size() == 1);

    /// Running one task at a time
    while (executable_task.size() != 0)
    {
        RAM* current_task;
        auto it = executable_task.begin();
        current_task = *it;
        current_task->set_comm(mcomm);

        /// Initialize all relations
        std::unordered_map<relation*, bool> scc_relation = current_task->get_RAM_relations();
        for (std::unordered_map<relation*, bool>::iterator it = scc_relation.begin() ; it != scc_relation.end(); ++it)
        {
            relation* rel = it->first;
            rel->initialize_relation_in_scc(it->second);
        }

        std::vector<u32> history;

        /// if case is for rules (acopy and copy) that requires only one iteration
        /// else case is for join rules

        int counter = 0;
        //std::cout << "-------------------B------------------" << std::endl;
        current_task->print_all_relation();
        if (current_task->get_iteration_count() == 1)
        {
            if (mcomm.get_local_rank() == 0)
                std::cout << "-------------------" << current_task->get_id() << " ITERATION " << counter << " ------------------" << std::endl;
            current_task->execute_in_batches(batch_size, history, intern_map);

        }
        else
        {
            u64 delta_in_scc = 0;
            do
            {
                if (mcomm.get_local_rank() == 0)
                    std::cout << "-------------------" << current_task->get_id() << " ITERATION " << counter << " ------------------" << std::endl;

                current_task->execute_in_batches(batch_size, history, intern_map);
                delta_in_scc = history[history.size()-2];
                counter++;


                //if (current_task->get_id() == 1)
                //    current_task->print_all_relation();
            }
            while (delta_in_scc != 0);
        }

        current_task->local_insert_in_full();
        //std::cout << "Number of iterations in scc " << current_task->get_id() << " " << counter << std::endl;

        //std::cout << "-------------------A------------------" << std::endl;
        //current_task->print_all_relation();

        /// marks executable_task as finished
        update_task_graph(executable_task);

        /// loads new runnable task
        executable_task = list_of_runnable_tasks(tasks, taskgraph1);
    }

    return true;
}


