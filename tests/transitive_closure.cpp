#include "../src/parallel_RA_inc.h"


int main(int argc, char **argv)
{
    mpi_comm mcomm;
    mcomm.create(argc, argv);

    relation* rel_path_2_1_2 = new relation(2, 2, 257, "D_path_2_1_2", FULL);
    relation* rel_edge_2_1_2 = new relation(2, 2, 256, "D_edge_2_1_2", FULL);
    relation* rel_path_2_1 = new relation(1, 2, 257, "D_path_2_1", FULL);
    relation* rel_edge_2_2 = new relation(1, 2, 256, "D_edge_2_2", FULL);

    RAM* scc5063 = new RAM(false, "copy");
    scc5063->add_relation(rel_edge_2_1_2);
    scc5063->add_relation(rel_path_2_1_2);
    scc5063->add_rule(new parallel_copy(rel_path_2_1_2, rel_edge_2_1_2, FULL, {0, 1, -1}));

    RAM* scc5064 = new RAM(true, "join");
    scc5064->add_relation(rel_edge_2_2);
    scc5064->add_relation(rel_path_2_1);
    scc5064->add_relation(rel_path_2_1_2);
    scc5064->add_rule(new parallel_acopy(rel_path_2_1, rel_path_2_1_2, DELTA, {0, 2, 1}));
    scc5064->add_rule(new parallel_join(rel_path_2_1_2, rel_path_2_1, DELTA, rel_edge_2_2, FULL, {-1, -1, 1, -1, 0}));

    RAM* scc5065 = new RAM(false, "acopy");
    scc5065->add_relation(rel_edge_2_2);
    scc5065->add_relation(rel_edge_2_1_2);
    scc5065->add_rule(new parallel_acopy(rel_edge_2_2, rel_edge_2_1_2, FULL, {2, 0, 1}));

    LIE* lie = new LIE();
    lie->add_relation(rel_path_2_1_2);
    lie->add_relation(rel_edge_2_1_2);
    lie->add_relation(rel_path_2_1);
    lie->add_relation(rel_edge_2_2);

    lie->add_scc(scc5063);
    lie->add_scc(scc5064);
    lie->add_scc(scc5065);
    lie->add_scc_dependance(scc5063, scc5064);
    lie->add_scc_dependance(scc5065, scc5064);

    lie->set_comm(mcomm);
    lie->set_batch_size(1);
    lie->execute();

    mcomm.destroy();
    return 0;
}
