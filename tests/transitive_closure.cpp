#include "../src/parallel_RA_inc.h"


int main(int argc, char **argv)
{
    mpi_comm mcomm;
    mcomm.create(argc, argv);    


    relation* rel_path_2_1_2 = new relation(2, 2, 257, "D_path_2_1_2", FULL);
    relation* rel_edge_2_1_2 = new relation(2, 2, 256, "D_edge_2_1_2", FULL);
    relation* rel_path_2_1 = new relation(1, 2, 257, "D_path_2_1", FULL);
    relation* rel_edge_2_2 = new relation(1, 2, 256, "D_edge_2_2", FULL);

    RAM* scc5026 = new RAM(false);
    scc5026->add_relation(rel_edge_2_1_2, false);
    scc5026->add_relation(rel_path_2_1_2, true);
    scc5026->add_rule(new parallel_copy(rel_path_2_1_2, rel_edge_2_1_2, FULL, {0, 1, -1}));

    RAM* scc5027 = new RAM(true);
    scc5027->add_relation(rel_edge_2_2, false);
    scc5027->add_relation(rel_path_2_1, true);
    scc5027->add_relation(rel_path_2_1_2, true);
    scc5027->add_rule(new parallel_acopy(rel_path_2_1, rel_path_2_1_2, DELTA, {0, 2, 1}));
    scc5027->add_rule(new parallel_join(rel_path_2_1_2, rel_path_2_1, DELTA, rel_edge_2_2, FULL, {-1, -1, 1, -1, 0}));

    RAM* scc5028 = new RAM(false);
    scc5028->add_relation(rel_edge_2_2, true);
    scc5028->add_relation(rel_edge_2_1_2, true);
    scc5028->add_rule(new parallel_acopy(rel_edge_2_2, rel_edge_2_1_2, DELTA, {2, 0, 1}));

    LIE* lie = new LIE();
    lie->add_relation(rel_path_2_1_2);
    lie->add_relation(rel_edge_2_1_2);
    lie->add_relation(rel_path_2_1);
    lie->add_relation(rel_edge_2_2);
    lie->add_scc(scc5026);
    lie->add_scc(scc5027);
    lie->add_scc(scc5028);
    lie->add_scc_dependance(scc5026, scc5027);
    lie->add_scc_dependance(scc5028, scc5027);


    lie->set_comm(mcomm);
    lie->set_batch_size(1);


    lie->execute();

    rel_path_2_1->print();

    mcomm.destroy();
    return 0;
}
