/*
 * Logical Inferencing Engine (LIE)
 * Copyright (c) Sidharth Kumar, et al, see License.md
 */



#include "../parallel_RA_inc.h"
//#include <experimental/filesystem>


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
            std::cout << curr_relation->get_debug_id() << ": {" << curr_relation->get_arity() << "}. (" << global_total_facts[i] << " total facts)" << std::endl;
        total_facts = total_facts + global_total_facts[i];
    }

    if (mcomm.get_local_rank() == 0)
        std::cout << "Total facts across all relations " << total_facts << std::endl << std::endl;
}



void LIE::write_checkpoint_dump(int loop_counter, std::vector<int> executed_scc_id, int scc_id)
{
    std::string dir_name;
    dir_name = output_dir + "/checkpoint-" + std::to_string(scc_id) + "-" + std::to_string(loop_counter);

    std::string scc_metadata;
    scc_metadata = dir_name + "/scc_metadata";
    if (mcomm.get_local_rank() == 0)
    {
        //mkdir(dir_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        FILE *fp;
        fp = fopen(scc_metadata.c_str(), "w");
        for (int i = 0; i < (int)executed_scc_id.size(); i++)
            fprintf (fp, "%d\n", executed_scc_id[i]);
        fclose(fp);
    }
    MPI_Barrier(mcomm.get_local_comm());

    for (u32 i = 0 ; i < lie_relation_count; i++)
    {
        //std::cout << "Name " << lie_relations[i]->get_filename() << std::endl;
        //if(lie_relations[i]->get_debug_id() == "rel_path_2_1_2")
        lie_relations[i]->parallel_IO(dir_name);
    }
}



void LIE::create_checkpoint_dump(int loop_counter, int scc_id)
{
    std::string dir_name;
    dir_name = output_dir + "/checkpoint-" + std::to_string(scc_id) + "-" + std::to_string(loop_counter);

    std::string scc_metadata;
    scc_metadata = dir_name + "/scc_metadata";
    if (mcomm.get_local_rank() == 0)
        mkdir(dir_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    MPI_Barrier(mcomm.get_local_comm());
}



bool LIE::execute ()
{
    /// Main : Execute : init : start
    mcomm.set_local_comm(MPI_COMM_WORLD);

#ifdef GOOGLE_MAP
    if (mcomm.get_local_rank() == 0)
        std::cout << "Using Google maps"  << std::endl;
#else
    if (mcomm.get_local_rank() == 0)
        std::cout << "Using Tom's Shmap" << std::endl;
#endif

    /// Initialize all relations
    for (u32 i = 0 ; i < lie_relation_count; i++)
    {
        lie_relations[i]->set_restart_flag(restart_flag);
        lie_relations[i]->set_share_io(share_io);
        lie_relations[i]->set_separate_io(separate_io);
        lie_relations[i]->set_offset_io(offset_io);
        lie_relations[i]->initialize_relation(mcomm);

#if DEBUG_OUTPUT
        //lie_relations[i]->print();
#endif
    }

    /// create output directory for checkpoint dumps
    if (enable_io == true && mcomm.get_local_rank() == 0)
        mkdir(output_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

#if DEBUG_OUTPUT
    if (mcomm.get_local_rank() == 0)
    {
        if (restart_flag == true)
        {
            if (separate_io == true)
                std::cout << "Read Local Data without MetaData" << std::endl;
            else if (offset_io == true)
                std::cout << "Read Global Data with Offset MetaData" << std::endl;
            else
                std::cout << "Read Global Data with Size MetaData" << std::endl;
        }
        else
        {
            if (separate_io == true)
                std::cout << "Write Local Data without MetaData" << std::endl;
            else if (offset_io == true)
                std::cout << "Write Global Data with Offset MetaData" << std::endl;
            else
                std::cout << "Write Global Data with Size MetaData" << std::endl;
        }
        std::cout << "----------------- Initialization Complete ---------------------" << std::endl << std::endl;
    }
#endif

    /// Executable task
    RAM* executable_task = one_runnable_tasks();



    std::vector<int> executed_scc_id;  /* the sccs that have been executed */
    /* Read scc metadate if it restarts from checkpoint */
    if (restart_flag == true)
    {
        int a;
        std::string scc_metadata;
        scc_metadata = restart_dir_name + "/scc_metadata";
        std::ifstream file(scc_metadata);
        if (file.is_open())
        {
            while (file >> a)
                executed_scc_id.push_back(a);
        }
        file.close();
    }

    all_to_all_buffer_size = new int**[lie_sccs_count];
    for (u32 i =0; i <lie_sccs_count; i++)
    {
        all_to_all_buffer_size[i] = new int*[MAX_LOOP_COUNT];
        for (u32 j =0; j <MAX_LOOP_COUNT; j++)
            all_to_all_buffer_size[i][j] = new int[mcomm.get_local_nprocs()];
    }

    compute_size1 = new int*[lie_sccs_count];
    for (u32 i =0; i <lie_sccs_count; i++)
    {
        compute_size1[i] = new int[MAX_LOOP_COUNT];
        memset(compute_size1[i], 0, MAX_LOOP_COUNT * sizeof(u32));
    }

    compute_size2 = new int*[lie_sccs_count];
    for (u32 i =0; i <lie_sccs_count; i++)
    {
        compute_size2[i] = new int[MAX_LOOP_COUNT];
        memset(compute_size2[i], 0, MAX_LOOP_COUNT * sizeof(u32));
    }

    double** running_time = new double*[lie_sccs_count];
    double** running_intra_bucket_comm = new double*[lie_sccs_count];
    double** running_buffer_allocate = new double*[lie_sccs_count];
    double** running_local_compute = new double*[lie_sccs_count];
    double** running_all_to_all = new double*[lie_sccs_count];
    double** running_buffer_free = new double*[lie_sccs_count];
    double** running_insert_newt = new double*[lie_sccs_count];
    double** running_insert_in_full = new double*[lie_sccs_count];
    double** running_fp = new double*[lie_sccs_count];

    for (u32 i =0; i <lie_sccs_count; i++)
    {
        running_time[i] = new double[MAX_LOOP_COUNT];
        running_intra_bucket_comm[i] = new double[MAX_LOOP_COUNT];
        running_buffer_allocate[i] = new double[MAX_LOOP_COUNT];
        running_local_compute[i] = new double[MAX_LOOP_COUNT];
        running_all_to_all[i] = new double[MAX_LOOP_COUNT];
        running_buffer_free[i] = new double[MAX_LOOP_COUNT];
        running_insert_newt[i] = new double[MAX_LOOP_COUNT];
        running_insert_in_full[i] = new double[MAX_LOOP_COUNT];
        running_fp[i] = new double[MAX_LOOP_COUNT];

        memset(running_time[i], 0, MAX_LOOP_COUNT * sizeof(double));
        memset(running_intra_bucket_comm[i], 0, MAX_LOOP_COUNT * sizeof(double));
        memset(running_buffer_allocate[i], 0, MAX_LOOP_COUNT * sizeof(double));
        memset(running_local_compute[i], 0, MAX_LOOP_COUNT * sizeof(double));
        memset(running_all_to_all[i], 0, MAX_LOOP_COUNT * sizeof(double));
        memset(running_buffer_free[i], 0, MAX_LOOP_COUNT * sizeof(double));
        memset(running_insert_newt[i], 0, MAX_LOOP_COUNT * sizeof(double));
        memset(running_insert_in_full[i], 0, MAX_LOOP_COUNT * sizeof(double));
        memset(running_fp[i], 0, MAX_LOOP_COUNT * sizeof(double));
    }



    /// Running one task at a time
    while (executable_task != NULL)
    {
        int loop_counter = 0;
        /* Skip the scc if it has been executed before */
        int scc_id = executable_task->get_id();
        std::vector<int>::iterator it;
        it = find (executed_scc_id.begin(), executed_scc_id.end(), scc_id);
        if (it != executed_scc_id.end())
        {
            update_task_graph(executable_task);
            executable_task = one_runnable_tasks();
            continue;
        }

        executable_task->set_comm(mcomm);

        /// Initialize all relations
        relation** scc_relation = executable_task->get_RAM_relations();
        bool* scc_relation_status = executable_task->get_RAM_relations_status();;
        u32 scc_relation_count = executable_task->get_ram_relation_count();
        if (restart_flag == false)
        {
            for (u32 i=0; i < scc_relation_count; i++)
            {
                if (scc_relation_status[i] == true)
                    scc_relation[i]->insert_full_in_delta();
            }
        }
        else
        {
            for (u32 i = 0 ; i < scc_relation_count; i++)
            {
                std::string delta_filename;
                delta_filename = restart_dir_name + "/" + scc_relation[i]->get_debug_id().c_str() + "_delta";

                if (separate_io == true)
                    delta_filename = delta_filename + "_" + std::to_string(mcomm.get_local_rank());

                scc_relation[i]->set_filename(delta_filename);
                scc_relation[i]->set_initailization_type(0);

                int is_access = access(delta_filename.c_str(), F_OK);
                int access_sum = 0;
                MPI_Allreduce(&is_access, &access_sum, 1, MPI_INT, MPI_SUM, mcomm.get_local_comm());

                if (access_sum + mcomm.get_local_nprocs() > 0)
                {
                    if (separate_io == true)
                        scc_relation[i]->load_data_from_separate_files();
                    else if (offset_io == true)
                        scc_relation[i]->load_data_from_file_with_offset();
                    else
                        scc_relation[i]->load_data_from_file();
                }
            }
        }

        std::vector<u32> history;

#if DEBUG_OUTPUT
        if (mcomm.get_local_rank() == 0)
            std::cout << "-------------------Executing SCC " << executable_task->get_id() << "------------------" << std::endl;
#endif

        /// if case is for rules (acopy and copy) that requires only one iteration
        /// else case is for join rules


        double writing_checkpoint_dump_time = 0;
        int checkpoint_dumps_num = 0;

        /// For SCCs that runs for only one iteration
        if (executable_task->get_iteration_count() == 1)
        {
            if (enable_data_io == true && loop_counter % cp_iteration == 0)
                create_checkpoint_dump(loop_counter, executable_task->get_id());

            //executable_task->execute_in_batches_with_all_to_all_threshold(batch_size, history, intern_map, &running_time, &running_intra_bucket_comm, &running_buffer_allocate, &running_local_compute, &running_all_to_all, &running_buffer_free, &running_insert_newt, &running_insert_in_full);
            executable_task->execute_in_batches(app_name, batch_size, history, intern_map, running_time, running_intra_bucket_comm, running_buffer_allocate, running_local_compute, running_all_to_all, running_buffer_free, running_insert_newt, running_insert_in_full, running_fp, loop_counter, executable_task->get_id(), output_dir, all_to_all_meta_data_dump, all_to_all_buffer_size, compute_size1, compute_size2);

            executed_scc_id.push_back(executable_task->get_id());
            if (enable_data_io == true && loop_counter % cp_iteration == 0)
            {
                double write_cp_start = MPI_Wtime();
                write_checkpoint_dump(loop_counter, executed_scc_id, executable_task->get_id());
                double write_cp_end = MPI_Wtime();
                writing_checkpoint_dump_time = (write_cp_end - write_cp_start);
                double max_write_cp_time = 0;
                MPI_Reduce(&writing_checkpoint_dump_time, &max_write_cp_time, 1, MPI_DOUBLE, MPI_MAX, 0, mcomm.get_comm());
                //if (mcomm.get_local_rank() == 0)
                //	std::cout << "Writing checkpoint dump " << checkpoint_dumps_num << " takes " << max_write_cp_time << "(s)" << std::endl;
                checkpoint_dumps_num++;
            }
            loop_counter++;
            iteration_count[executable_task->get_id()] = loop_counter;


#if DEBUG_OUTPUT
            //for (u32 i = 0 ; i < scc_relation_count; i++)
            //    scc_relation[i]->print();
            print_all_relation_size();
#endif
        }
        /// For SCCs that runs till fixed point is reached
        else
        {
            if (mcomm.get_rank() == 0)
                std::cout << "Iteration#\tBuffer_creation_time\tComputation_time\tAll_to_all_time\tBuffer_free_time\tInsert_in_newt_time\tIntra_comm_time\tInsert_in_full_time\tTotal_time" << std::endl;
            u64 delta_in_scc = 0;
            do
            {
                if (enable_data_io == true && loop_counter % cp_iteration == 0)
                    create_checkpoint_dump(loop_counter, executable_task->get_id());

                //executable_task->execute_in_batches_with_all_to_all_threshold(batch_size, history, intern_map, &running_time, &running_intra_bucket_comm, &running_buffer_allocate, &running_local_compute, &running_all_to_all, &running_buffer_free, &running_insert_newt, &running_insert_in_full);
                executable_task->execute_in_batches(app_name, batch_size, history, intern_map, running_time, running_intra_bucket_comm, running_buffer_allocate, running_local_compute, running_all_to_all, running_buffer_free, running_insert_newt, running_insert_in_full, running_fp, loop_counter, executable_task->get_id(), output_dir, all_to_all_meta_data_dump, all_to_all_buffer_size, compute_size1, compute_size2);
                //executable_task->print_all_relation();

                delta_in_scc = history[history.size()-2];
                if (delta_in_scc == 0)
                    executed_scc_id.push_back(executable_task->get_id());
                if (enable_data_io == true && loop_counter % cp_iteration == 0)
                {
                    double write_cp_start = MPI_Wtime();
                    write_checkpoint_dump(loop_counter, executed_scc_id, executable_task->get_id());
                    double write_cp_end = MPI_Wtime();
                    writing_checkpoint_dump_time = (write_cp_end - write_cp_start);
                    double max_write_cp_time = 0;
                    MPI_Reduce(&writing_checkpoint_dump_time, &max_write_cp_time, 1, MPI_DOUBLE, MPI_MAX, 0, mcomm.get_comm());
                    //if (mcomm.get_local_rank() == 0)
                    //    std::cout << "Writing checkpoint dump " << checkpoint_dumps_num << " takes " << max_write_cp_time << "(s)" << std::endl;
                    checkpoint_dumps_num++;
                }
                loop_counter++;
                iteration_count[executable_task->get_id()] = loop_counter;


#if DEBUG_OUTPUT
                //for (u32 i = 0 ; i < scc_relation_count; i++)
                //    scc_relation[i]->print();
                //print_all_relation_size();
#endif
            }
            while (delta_in_scc != 0);
        }

        set_executed_scc_id(executed_scc_id);


        set_loop_counter(loop_counter);


        executable_task->insert_delta_in_full();

        /// marks executable_task as finished
        update_task_graph(executable_task);



        /// loads new runnable task
        executable_task = one_runnable_tasks();
    }

#if 1
    if (all_to_all_meta_data_dump == true)
    {
        double total_time = 0;
        double max_total_time = 0;
        for (u32 l=0; l < lie_sccs_count; l++)
            for (u32 il=0; il < iteration_count[l]; il++)
                total_time = total_time + running_time[l][il];

        MPI_Allreduce(&total_time, &max_total_time, 1, MPI_DOUBLE, MPI_MAX, mcomm.get_local_comm());
        if (total_time == max_total_time)
        {
            std::string time_file;
            time_file = output_dir + "/MAX_TIME_" + std::to_string(mcomm.get_rank());
            FILE *fpm;
            fpm = fopen(time_file.c_str(), "w");

            fprintf(fpm, "TaskID\tIteration#\tBuffer_creation_time\tComputation_time\tAll_to_all_time\tBuffer_free_time\tInsert_in_newt_time\tIntra_comm_time\tInsert_in_full_time\tFixed_point_check\tTotal_time\n");

            for (u32 l=0; l < lie_sccs_count; l++)
            {
                for (u32 il=0; il < iteration_count[l]; il++)
                {
                    fprintf(fpm, "%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", l, il, running_buffer_allocate[l][il], running_local_compute[l][il], running_all_to_all[l][il], running_buffer_free[l][il], running_insert_newt[l][il], running_intra_bucket_comm[l][il], running_insert_in_full[l][il], running_fp[l][il], running_time[l][il]);
                }
            }
            fclose(fpm);
        }

        if (mcomm.get_local_rank() == 0)
        {
            std::string nprocs_file;
            nprocs_file = output_dir + "/metadata";
            FILE *fpn;
            fpn = fopen(nprocs_file.c_str(), "w");
            fprintf(fpn, "(nprocs)\n%d\n", mcomm.get_local_nprocs());
            fprintf(fpn, "(task count)\n%d\n", lie_sccs_count);
            for (u32 l=0; l < lie_sccs_count; l++)
               fprintf(fpn, "(task %d)\n%d\n", l, iteration_count[l]);
            fclose(fpn);
        }

        for (u32 l=0; l < lie_sccs_count; l++)
        {
            for (u32 il=0; il < iteration_count[l]; il++)
            {
                if (enable_io == true)
                    create_checkpoint_dump(il, l);

                std::string dir_name;
                dir_name = output_dir + "/checkpoint-" + std::to_string(l) + "-" + std::to_string(il);
                std::string scc_metadata_a2a;
                scc_metadata_a2a = dir_name + "/a2a_" + std::to_string(mcomm.get_local_rank());

                FILE *fp1;
                fp1 = fopen(scc_metadata_a2a.c_str(), "w");
                for (int i=0; i< mcomm.get_local_nprocs(); i++)
                    fprintf (fp1, "%d\n", all_to_all_buffer_size[l][il][i]);
                fclose(fp1);


                std::string scc_metadata_lc;
                scc_metadata_lc = dir_name + "/lc";

                FILE *fp2;
                fp2 = fopen(scc_metadata_lc.c_str(), "w");
                fprintf (fp2, "%d\n", compute_size1[l][il]);
                fprintf (fp2, "%d\n", compute_size2[l][il]);
                fprintf (fp2, "%f\n", running_buffer_allocate[l][il]);
                fprintf (fp2, "%f\n", running_local_compute[l][il]);
                fprintf (fp2, "%f\n", running_all_to_all[l][il]);
                fprintf (fp2, "%f\n", running_buffer_free[l][il]);
                fprintf (fp2, "%f\n", running_insert_newt[l][il]);
                fprintf (fp2, "%f\n", running_intra_bucket_comm[l][il]);
                fprintf (fp2, "%f\n", running_insert_in_full[l][il]);
                fprintf (fp2, "%f\n", running_fp[l][il]);
                fprintf (fp2, "%f\n", running_time[l][il]);
                fclose(fp2);
            }
        }
    }


    for (u32 i =0; i <lie_sccs_count; i++)
    {
        for (u32 j =0; j <MAX_LOOP_COUNT; j++)
            delete[] all_to_all_buffer_size[i][j];
        delete[] all_to_all_buffer_size[i];
    }
    delete[] all_to_all_buffer_size;

    for (u32 i =0; i <lie_sccs_count; i++)
        delete[] compute_size1[i];
    delete[] compute_size1;

    for (u32 i =0; i <lie_sccs_count; i++)
        delete[] compute_size2[i];
    delete[] compute_size2;

    for (u32 i =0; i <lie_sccs_count; i++)
    {
        delete[] running_time[i];
        delete[] running_intra_bucket_comm[i];
        delete[] running_buffer_allocate[i];
        delete[] running_local_compute[i];
        delete[] running_all_to_all[i];
        delete[] running_buffer_free[i];
        delete[] running_insert_newt[i];
        delete[] running_insert_in_full[i];
        delete[] running_fp[i];
    }

    delete[] running_time;
    delete[] running_intra_bucket_comm;
    delete[] running_buffer_allocate;
    delete[] running_local_compute;
    delete[] running_all_to_all;
    delete[] running_buffer_free;
    delete[] running_insert_newt;
    delete[] running_insert_in_full;
    delete[] running_fp;
#endif





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
