#include "../src/parallel_RA_inc.h"


int main(int argc, char **argv)
{
    mpi_comm mcomm;
    mcomm.create(argc, argv);

    if (argc != 4 && argc != 5 && argc != 2)
    {
        printf("argc %d\n", argc);
        std::cout << "1. Usage: mpirun -n <process cout> ./TC <input data file> <output_dir> <iterations for checkpoint>\n"
                     "2. Restart: mpirun -n <process cout> ./TC --restart <checkpoint dump directory> <output_dir> <iterations for checkpoint>" << std::endl;
        MPI_Abort(mcomm.get_comm(), -1);
    }

    char rel_path_212[1024];
    char rel_edge_212[1024];
    char rel_path_21[1024];
    char rel_edge_22[1024];

    bool restart_flag = false;
    char* dir_name;
    char* output_dir;
    int cp_iteration;
    if (argc == 5 && strcmp(argv[1], "--restart") == 0)
    {
        restart_flag = true;
        dir_name = argv[2];
        output_dir = argv[3];
        cp_iteration = atoi(argv[4]);
        if ( dir_name[strlen(dir_name) - 1] == '/' ) dir_name[strlen(dir_name) - 1] = '\0';
        sprintf(rel_path_212, "%s/%s_full", dir_name, "T");
        sprintf(rel_edge_212, "%s/%s_full", dir_name, "G");
        sprintf(rel_path_21, "%s/%s_full", dir_name, "rel_path_2_1");
        sprintf(rel_edge_22, "%s/%s_full", dir_name, "rel_edge_2_2");
    }
    else if (argc == 4)
    {
        output_dir = argv[2];
        cp_iteration = atoi(argv[3]);
        sprintf(rel_path_212, "%s", "../data/g5955/path_2_1_2");
        sprintf(rel_edge_212, "%s", argv[1]);
    }
    else if (argc == 2)
    {
        sprintf(rel_path_212, "%s", "../data/g5955/path_2_1_2");
        sprintf(rel_edge_212, "%s", argv[1]);
    }
    else
    {
        std::cout << "1. Usage: mpirun -n <process cout> ./TC <input data file>\n"
                     "2. Restart: mpirun -n <process cout> ./TC --restart <checkpoint dump directory>" << std::endl;
    }

    relation* T = new relation(1, true, 2, 257, "T", rel_path_212, FULL);
    relation* G = new relation(1, true, 2, 256, "G", rel_edge_212, FULL);



    RAM* scc13237 = new RAM(true, 1);
    scc13237->add_relation(G, false);
    scc13237->add_relation(T, true);
    scc13237->add_rule(new parallel_join(T, T, DELTA, G, FULL, {3, 1}));


    RAM* scc13239 = new RAM(false, 2);
    scc13239->add_relation(G, false);
    scc13239->add_relation(T, true);
    scc13239->add_rule(new parallel_copy(T, G, FULL, {1, 0}));

    LIE* lie = new LIE();
    lie->add_relation(T);
    lie->add_relation(G);
    lie->add_scc(scc13237);
    lie->add_scc(scc13239);
    lie->add_scc_dependance(scc13239, scc13237);

    if (restart_flag == true)
        lie->set_restart_dir_name(dir_name);


    if (argc != 2)
    {
        lie->enable_all_to_all_dump();
        lie->enable_separate_io();

        lie->set_cp_iteration(cp_iteration);
        lie->set_output_dir(output_dir);
        lie->set_restart_flag(restart_flag); // set restart flag
        lie->enable_IO();
    }

    lie->set_name(argv[0]);
    lie->set_comm(mcomm);
    lie->set_batch_size(1);
    lie->execute();
    lie->print_all_relation_size();

#if 0
    if (argc != 2)
    {
        double write_cp_start = MPI_Wtime();
        int loop_count = lie->get_loop_counter();
        std::vector<int> executed_scc_id = lie->get_executed_scc_id();
        lie->write_checkpoint_dump(loop_count, executed_scc_id);
        double write_cp_end = MPI_Wtime();
        double writing_checkpoint_dump_time = (write_cp_end - write_cp_start);
        double max_write_cp_time = 0;
        MPI_Reduce(&writing_checkpoint_dump_time, &max_write_cp_time, 1, MPI_DOUBLE, MPI_MAX, 0, mcomm.get_comm());
        if (mcomm.get_rank() == 0)
            std::cout << "Writing last checkpoint dump takes " << max_write_cp_time << "(s)" << std::endl;
    }
#endif

    delete lie;

    mcomm.destroy();
    return 0;
}
