/*
 * Logical Inferencing Engine (LIE)
 * Copyright (c) Sidharth Kumar, et al, see License.md
 */



#include "../parallel_RA_inc.h"


void LIE::add_relation(relation* rel)
{
    lie_relations[lie_relation_count] = rel;
    lie_relation_count++;
}

void LIE::add_scc(RAM* ra)
{
    lie_sccs[lie_sccs_count] = ra;
    lie_sccs_count++;
}


/// This function currentl only returns one runnable task
/// size of returned list is always going to be 1
RAM* LIE::one_runnable_tasks()
{
    u32 counter = 0;
    for (u32 i=0; i < lie_sccs_count; i++)
    {
        if (lie_sccs[i] == NULL)
        {
            counter++;
            continue;
        }
        if (counter == lie_sccs_count)
            return NULL;

        bool break_loop = false;
        for (auto it = taskgraph.begin(); it != taskgraph.end(); it++)
        {
            std::set<RAM*> it2 = it->second;
            for (auto dit2 = it2.begin(); dit2 != it2.end(); dit2++)
            {
                if (lie_sccs[i] == *dit2)
                {
                    break_loop=true;
                    break;
                }
            }
            if (break_loop==true)
                break;
        }
        if (break_loop==false)
            return lie_sccs[i];
    }

    return NULL;
}



void LIE::update_task_graph(RAM* executable_task)
{
    for (u32 i=0; i < lie_sccs_count; i++)
    {
        if (lie_sccs[i] == executable_task)
        {
            taskgraph.erase(lie_sccs[i]);
            delete lie_sccs[i];
            lie_sccs[i] = NULL;
        }
    }
}



void LIE::add_scc_dependance (RAM* src_task, RAM* destination_task)
{
    auto it = taskgraph.find(src_task);
    if( it != taskgraph.end() )
    {
        auto it2 = (it->second).find(destination_task);
        if( it2 == (it->second).end() )
        {
            (it->second).insert(destination_task);
            taskgraph[src_task] = it->second;
        }
    }
    else
    {
        std::set<RAM*> k;
        k.insert(destination_task);
        taskgraph.insert(std::make_pair(src_task, k));
    }

    std::set<RAM*> temp = taskgraph[src_task];
    temp.insert(destination_task);

    return;
}


#if 0
void LIE::print_all_relation()
{
    u64 total_facts=0;
    for (u32 i = 0 ; i < lie_relation_count; i++)
    {
        relation* curr_relation = lie_relations[i];
        u64 local_facts = curr_relation->get_full_element_count();
        u64 global_total_facts = 0;
        MPI_Allreduce(&local_facts, &global_total_facts, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mcomm.get_local_comm());

        if (mcomm.get_local_rank() == 0)
            std::cout << curr_relation->get_debug_id() << ": {" << curr_relation->get_arity() << "}. (" << global_total_facts << " total facts)" << std::endl;

        total_facts = total_facts + global_total_facts;
    }
    if (mcomm.get_local_rank() == 0)
        std::cout << "Total facts across all relations " << total_facts << std::endl;
}
#endif



void LIE::print_all_relation_size()
{
    u64 total_facts=0;
    u64 local_facts[lie_relation_count];
    for (u32 i = 0 ; i < lie_relation_count; i++)
    {
        relation* curr_relation = lie_relations[i];
        local_facts[i] = curr_relation->get_full_element_count();
    }

    u64 global_total_facts[lie_relation_count];
    MPI_Allreduce(local_facts, global_total_facts, lie_relation_count, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mcomm.get_local_comm());
    for (u32 i = 0 ; i < lie_relation_count; i++)
    {
        relation* curr_relation = lie_relations[i];
        if (mcomm.get_local_rank() == 0)
        {
            std::cout << curr_relation->get_debug_id() << ": {" << curr_relation->get_arity() << "}. (" << global_total_facts[i] << " total facts)" << std::endl;
            //debug_buffer.push_back(curr_relation->get_debug_id() + ": {" + std::to_string(curr_relation->get_arity()) + "}. (" + std::to_string(global_total_facts[i]) + " total facts)\n");
        }
        total_facts = total_facts + global_total_facts[i];
    }

    if (mcomm.get_local_rank() == 0)
    {
        //debug_buffer.push_back("Total facts across all relations " + std::to_string(total_facts) + "\n\n");
        std::cout << "Total facts across all relations " << total_facts << std::endl << std::endl;
    }
}



bool LIE::execute ()
{
    /// Main : Execute : init : start
    mcomm.set_local_comm(MPI_COMM_WORLD);
    debug_buffer.reserve(400000);

    /// Initialize all relations
    for (u32 i = 0 ; i < lie_relation_count; i++)
    {
        lie_relations[i]->initialize_relation(mcomm);

#if DEBUG_OUTPUT
        //lie_relations[i]->print();
#endif
    }


#if DEBUG_OUTPUT
    if (mcomm.get_local_rank() == 0)
    {
        //std::cout << "----------------- Initialization Complete ---------------------" << std::endl << std::endl;
        debug_buffer.push_back("----------------- Initialization Complete ---------------------\n\n");
    }
#endif


    /// Executable task
    RAM* executable_task = one_runnable_tasks();

    double running_time=0;
    double running_intra_bucket_comm=0;
    double running_buffer_allocate=0;
    double running_local_compute=0;
    double running_all_to_all=0;
    double running_buffer_free=0;
    double running_insert_newt=0;
    double running_insert_in_full=0;
    double running_fp = 0;
    int loop_counter = 0;

    /// Running one task at a time
    while (executable_task != NULL)
    {
        executable_task->set_comm(mcomm);

        /// Initialize all relations
        relation** scc_relation = executable_task->get_RAM_relations();
        bool* scc_relation_status = executable_task->get_RAM_relations_status();;
        u32 scc_relation_count = executable_task->get_ram_relation_count();
        for (u32 i=0; i < scc_relation_count; i++)
        {
            if (scc_relation_status[i] == true)
                scc_relation[i]->insert_full_in_delta();
        }

        std::vector<u32> history;

#if DEBUG_OUTPUT
        if (mcomm.get_local_rank() == 0)
        {
            //std::cout << "-------------------Executing SCC " << executable_task->get_id() << "------------------" << std::endl;
            debug_buffer.push_back("-------------------Executing SCC " + std::to_string(executable_task->get_id()) + " ------------------\n");
        }
#endif

        /// if case is for rules (acopy and copy) that requires only one iteration
        /// else case is for join rules
        /// For SCCs that runs for only one iteration
        if (executable_task->get_iteration_count() == 1)
        {
            executable_task->execute_in_batches(batch_size, history, intern_map, &running_time, &running_intra_bucket_comm, &running_buffer_allocate, &running_local_compute, &running_all_to_all, &running_buffer_free, &running_insert_newt, &running_insert_in_full, &running_fp, debug_buffer, debug_file_name);
            loop_counter++;


#if DEBUG_OUTPUT
            //for (u32 i = 0 ; i < scc_relation_count; i++)
            //    scc_relation[i]->print();
            print_all_relation_size();
#endif
        }

        /// For SCCs that runs till fixed point is reached
        else
        {
            u64 delta_in_scc = 0;
            do
            {
                executable_task->execute_in_batches(batch_size, history, intern_map, &running_time, &running_intra_bucket_comm, &running_buffer_allocate, &running_local_compute, &running_all_to_all, &running_buffer_free, &running_insert_newt, &running_insert_in_full, &running_fp, debug_buffer, debug_file_name);
                loop_counter++;
                delta_in_scc = history[history.size()-2];


                if (enable_io == true)
                {
                    if (loop_counter % 10 == 0)
                    {
                        char dir_name[1024];
                        sprintf(dir_name, "output/checkpoin-%d", loop_counter);
                        if (mcomm.get_local_rank() == 0)
                        {
                            mkdir("output", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                            mkdir(dir_name, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                        }
                        MPI_Barrier(mcomm.get_local_comm());

                        for (u32 i = 0 ; i < lie_relation_count; i++)
                            lie_relations[i]->parallel_IO(dir_name);
                    }
                }

#if DEBUG_OUTPUT
                //for (u32 i = 0 ; i < scc_relation_count; i++)
                //    scc_relation[i]->print();
                print_all_relation_size();
#endif
            }
            while (delta_in_scc != 0);
        }


        executable_task->insert_delta_in_full();

        /// marks executable_task as finished
        update_task_graph(executable_task);

        /// loads new runnable task
        executable_task = one_runnable_tasks();
    }

    //std::ofstream myfile;
    //myfile.open (debug_file_name);
    //for (std::string n : debug_buffer)
    //    myfile << n;
    //myfile.close();

    return true;
}


LIE::~LIE ()
{
    for (u32 i = 0 ; i < lie_relation_count; i++)
    {
        lie_relations[i]->finalize_relation();
        delete (lie_relations[i]);
    }
}
