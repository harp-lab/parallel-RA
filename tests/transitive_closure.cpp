#include "../src/parallel_RA_inc.h"


int main(int argc, char **argv)
{
    mpi_comm mcomm;
    mcomm.create(argc, argv);    

    if (argc != 3 && argc != 4)
    {
    	printf("argc %d\n", argc);
        std::cout << "1. Usage: mpirun -n <process cout> ./TC <input data file> <output_dir>\n"
        		"2. Restart: mpirun -n <process cout> ./TC --restart <checkpoint dump directory> <output_dir>" << std::endl;
        MPI_Abort(mcomm.get_comm(), -1);
    }

    char rel_path_212[1024];
    char rel_edge_212[1024];
    char rel_path_21[1024];
    char rel_edge_22[1024];

    bool restart_flag = false;
    char* dir_name;
    char* output_dir;
    if (argc == 4 && strcmp(argv[1], "--restart") == 0)
    {
    	restart_flag = true;
    	dir_name = argv[2];
    	output_dir = argv[3];
    	if ( dir_name[strlen(dir_name) - 1] == '/' ) dir_name[strlen(dir_name) - 1] = '\0';
    	sprintf(rel_path_212, "%s/%s_full", dir_name, "rel_path_2_1_2");
    	sprintf(rel_edge_212, "%s/%s_full", dir_name, "rel_edge_2_1_2");
    	sprintf(rel_path_21, "%s/%s_full", dir_name, "rel_path_2_1");
    	sprintf(rel_edge_22, "%s/%s_full", dir_name, "rel_edge_2_2");
    }
    else if (argc == 3)
    {
    	output_dir = argv[2];
    	sprintf(rel_path_212, "%s", "../data/g5955/path_2_1_2");
		sprintf(rel_edge_212, "%s", argv[1]);
		sprintf(rel_path_21, "%s", "../data/g5955/path_2_1");
		sprintf(rel_edge_22, "%s", "../data/g5955/edge_2_2");
    }
    else
    {
    	std::cout << "1. Usage: mpirun -n <process cout> ./TC <input data file>\n"
    	        "2. Restart: mpirun -n <process cout> ./TC --restart <checkpoint dump directory>" << std::endl;
    }

    relation* rel_path_2_1_2 = new relation(2, true, 2, 257, "rel_path_2_1_2", rel_path_212, FULL);
    relation* rel_edge_2_1_2 = new relation(2, true, 2, 256, "rel_edge_2_1_2", rel_edge_212, FULL);
    relation* rel_path_2_1 = new relation(1, false, 2, 257, "rel_path_2_1", rel_path_21, FULL);
    relation* rel_edge_2_2 = new relation(1, false, 2, 256, "rel_edge_2_2", rel_edge_22, FULL);

    RAM* scc13237 = new RAM(true, 1);
    scc13237->add_relation(rel_edge_2_2, false);
    scc13237->add_relation(rel_path_2_1, true);
    scc13237->add_relation(rel_path_2_1_2, true);
    scc13237->add_rule(new parallel_acopy(rel_path_2_1, rel_path_2_1_2, DELTA, {0, 2, 1}));
    scc13237->add_rule(new parallel_join(rel_path_2_1_2, rel_path_2_1, DELTA, rel_edge_2_2, FULL, {4, 2}));

    RAM* scc13238 = new RAM(false, 3);
    scc13238->add_relation(rel_edge_2_2, true);
    scc13238->add_relation(rel_edge_2_1_2, true);
    scc13238->add_rule(new parallel_acopy(rel_edge_2_2, rel_edge_2_1_2, DELTA, {1, 2, 0}));

    RAM* scc13239 = new RAM(false, 2);
    scc13239->add_relation(rel_edge_2_1_2, false);
    scc13239->add_relation(rel_path_2_1_2, true);
    scc13239->add_rule(new parallel_copy(rel_path_2_1_2, rel_edge_2_1_2, FULL, {0, 1}));

    LIE* lie = new LIE();
    lie->add_relation(rel_path_2_1_2);
    lie->add_relation(rel_edge_2_1_2);
    lie->add_relation(rel_path_2_1);
    lie->add_relation(rel_edge_2_2);
	lie->add_scc(scc13237);
	lie->add_scc(scc13238);
	lie->add_scc(scc13239);
	lie->add_scc_dependance(scc13238, scc13237);
	lie->add_scc_dependance(scc13239, scc13237);

	if (restart_flag == true)
		lie->set_restart_dir_name(dir_name);

//	lie->enable_offset_io();
//	lie->enable_separate_io();
	lie->enable_share_io();

	lie->set_output_dir(output_dir);
    lie->set_restart_flag(restart_flag); // set restart flag
    lie->enable_IO();
    lie->set_comm(mcomm);
    lie->set_batch_size(1);
    lie->execute();
    lie->print_all_relation_size();

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

    delete lie;

    mcomm.destroy();
    return 0;
}
