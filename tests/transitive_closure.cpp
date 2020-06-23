#include "../src/parallel_RA_inc.h"

// Incorporate load balancer
// fix all arity issues
// set get function issues
// logging stuff
// join and copy two time access issue
// untangle task level parallelism
// cleanup task level parallelism
// cleanup cumulative all to all
// document the code
// intern

int main(int argc, char **argv)
{
    mpi_comm mcomm;
    mcomm.create(argc, argv);

    //u32 batch_size = atoi(argv[1]);
    //double batch_time = atof(argv[2]);
    //float threshold = atof(argv[3]);
    //int mode = atoi(argv[4]);
    //u32 task_count = atoi(argv[5]);

    relation* rel_path_2_1_2 = new relation(2, 2, 257, "/var/tmp/g4831/path_2_1_2", FULL);
    relation* rel_edge_2_ = new relation(0, 2, 256, "/var/tmp/g4831/edge_2_", FULL);
    relation* rel_edge_2_1_2 = new relation(2, 2, 256, "/var/tmp/g4831/edge_2_1_2", FULL);
    relation* rel_path_2_1 = new relation(1, 2, 257, "/var/tmp/g4831/path_2_1", FULL);
    relation* rel_edge_2_2 = new relation(1, 2, 256, "/var/tmp/g4831/edge_2_2", FULL);

    RAM* scc4832 = new RAM();
    scc4832->add_relation(rel_path_2_1_2);
    scc4832->add_relation(rel_edge_2_);
    scc4832->add_relation(rel_edge_2_1_2);
    scc4832->add_relation(rel_path_2_1);
    scc4832->add_relation(rel_edge_2_2);
    scc4832->add_rule(new parallel_copy(rel_path_2_1_2, rel_edge_2_, FULL, {-1, 0, 1}));

    RAM* scc4833 = new RAM();
    scc4833->add_relation(rel_path_2_1_2);
    scc4833->add_relation(rel_edge_2_);
    scc4833->add_relation(rel_edge_2_1_2);
    scc4833->add_relation(rel_path_2_1);
    scc4833->add_relation(rel_edge_2_2);
    scc4833->add_rule(new parallel_acopy(rel_edge_2_2, rel_edge_2_1_2, DELTA, {2, 0, 1}));

    RAM* scc4834 = new RAM();
    scc4834->add_relation(rel_path_2_1_2);
    scc4834->add_relation(rel_edge_2_);
    scc4834->add_relation(rel_edge_2_1_2);
    scc4834->add_relation(rel_path_2_1);
    scc4834->add_relation(rel_edge_2_2);
    scc4834->add_rule(new parallel_acopy(rel_edge_2_, rel_edge_2_1_2, DELTA, {1, 2, 0}));

    RAM* scc4835 = new RAM();
    scc4835->add_relation(rel_path_2_1_2);
    scc4835->add_relation(rel_edge_2_);
    scc4835->add_relation(rel_edge_2_1_2);
    scc4835->add_relation(rel_path_2_1);
    scc4835->add_relation(rel_edge_2_2);
    scc4835->add_rule(new parallel_acopy(rel_path_2_1, rel_path_2_1_2, DELTA, {0, 2, 1}));
    scc4835->add_rule(new parallel_join(rel_path_2_1_2, rel_path_2_1, DELTA, rel_edge_2_2, FULL, {-1, -1, 1, -1, 0}, 1));

    LIE* lie = new LIE();
    lie->add_scc(scc4832);
    lie->add_scc(scc4833);
    lie->add_scc(scc4834);
    lie->add_scc(scc4835);
    lie->add_scc_dependance(scc4835, scc4832);
    lie->add_scc_dependance(scc4835, scc4833);
    lie->add_scc_dependance(scc4832, scc4834);
    lie->set_comm(mcomm);
    lie->set_batch_size(10);
    lie->execute();


    mcomm.destroy();
    return 0;
}
