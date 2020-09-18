#include "../src/parallel_RA_inc.h"


int main(int argc, char **argv)
{
    mpi_comm mcomm;
    mcomm.create(argc, argv);



relation* rel_4cl_4_4_3_2 = new relation(3, false, 4, 259, "rel_4cl_4_4_3_2", "../data/5clique//4cl_4_17", FULL);
relation* rel_2cl_2_2_1 = new relation(2, true, 2, 256, "rel_2cl_2_2_1", "../data/5clique//2cl_2_16", FULL);
relation* rel_inner_replacement21_4_1_2_3_4 = new relation(4, true, 4, 257, "rel_inner_replacement21_4_1_2_3_4", "../data/5clique//inner-replacement21_4_15", FULL);
relation* rel_4cl_4_3_2_1 = new relation(3, false, 4, 259, "rel_4cl_4_3_2_1", "../data/5clique//4cl_4_14", FULL);
relation* rel_4cl_4_1_2_3_4 = new relation(4, true, 4, 259, "rel_4cl_4_1_2_3_4", "../data/5clique//4cl_4_13", FULL);
relation* rel_5cl_5_1_2_3_4_5 = new relation(5, true, 5, 264, "rel_5cl_5_1_2_3_4_5", "../data/5clique//5cl_5_12", FULL);
relation* rel_2cl_2_1 = new relation(1, false, 2, 256, "rel_2cl_2_1", "../data/5clique//2cl_2_11", FULL);
relation* rel_inner_replacement21_4_4_1 = new relation(2, false, 4, 257, "rel_inner_replacement21_4_4_1", "../data/5clique//inner-replacement21_4_10", FULL);
relation* rel_2cl_2_2 = new relation(1, false, 2, 256, "rel_2cl_2_2", "../data/5clique//2cl_2_9", FULL);
relation* rel_graph_2_2_1 = new relation(2, true, 2, 258, "rel_graph_2_2_1", "../data/5clique//graph_2_8", FULL);
relation* rel_graph_2_1_2 = new relation(2, true, 2, 258, "rel_graph_2_1_2", "../data/5clique//graph_2_7", FULL);
relation* rel_inter_body26_3_1_2_3 = new relation(3, true, 3, 263, "rel_inter_body26_3_1_2_3", "../data/5clique//inter-body26_3_6", FULL);
relation* rel_3cl_3_1_2_3 = new relation(3, true, 3, 262, "rel_3cl_3_1_2_3", "../data/5clique//3cl_3_5", FULL);
relation* rel_inner_replacement7_5_1_2_3_4_5 = new relation(5, true, 5, 261, "rel_inner_replacement7_5_1_2_3_4_5", "../data/5clique//inner-replacement7_5_4", FULL);
relation* rel_3cl_3_2_1 = new relation(2, false, 3, 262, "rel_3cl_3_2_1", "../data/5clique//3cl_3_3", FULL);
relation* rel_inner_replacement7_5_5_1 = new relation(2, false, 5, 261, "rel_inner_replacement7_5_5_1", "../data/5clique//inner-replacement7_5_2", FULL);
relation* rel_inter_body26_3_1_3 = new relation(2, false, 3, 263, "rel_inter_body26_3_1_3", "../data/5clique//inter-body26_3_1", FULL);
relation* rel_3cl_3_3_2 = new relation(2, false, 3, 262, "rel_3cl_3_3_2", "../data/5clique//3cl_3_0", FULL);

RAM* scc4696 = new RAM(false, 1);
scc4696->add_relation(rel_graph_2_2_1, false);
scc4696->add_relation(rel_2cl_2_2_1, true);
scc4696->add_rule(new parallel_copy_filter(rel_2cl_2_2_1, rel_graph_2_2_1, FULL, {0, 1}, [](const u64* const data){ return ((s32)(data[0]&0xffffffff) < (s32)(data[1]&0xffffffff)); }));

RAM* scc4697 = new RAM(false, 5);
scc4697->add_relation(rel_3cl_3_3_2, true);
scc4697->add_relation(rel_3cl_3_1_2_3, true);
scc4697->add_rule(new parallel_acopy(rel_3cl_3_3_2, rel_3cl_3_1_2_3, DELTA, {2, 1, 3, 0}));

RAM* scc4698 = new RAM(false, 9);
scc4698->add_relation(rel_inner_replacement21_4_4_1, true);
scc4698->add_relation(rel_inner_replacement21_4_1_2_3_4, true);
scc4698->add_rule(new parallel_acopy(rel_inner_replacement21_4_4_1, rel_inner_replacement21_4_1_2_3_4, DELTA, {3, 0, 4, 1, 2}));

RAM* scc4699 = new RAM(false, 13);
scc4699->add_relation(rel_3cl_3_2_1, true);
scc4699->add_relation(rel_3cl_3_1_2_3, true);
scc4699->add_rule(new parallel_acopy(rel_3cl_3_2_1, rel_3cl_3_1_2_3, DELTA, {1, 0, 3, 2}));

RAM* scc4700 = new RAM(false, 18);
scc4700->add_relation(rel_inter_body26_3_1_2_3, true);
scc4700->add_relation(rel_2cl_2_2, false);
scc4700->add_relation(rel_2cl_2_1, false);
scc4700->add_rule(new parallel_join(rel_inter_body26_3_1_2_3, rel_2cl_2_2, FULL, rel_2cl_2_1, FULL, {4, 0, 2}));

RAM* scc4701 = new RAM(false, 3);
scc4701->add_relation(rel_inner_replacement7_5_5_1, false);
scc4701->add_relation(rel_5cl_5_1_2_3_4_5, true);
scc4701->add_relation(rel_2cl_2_2_1, false);
scc4701->add_rule(new parallel_join(rel_5cl_5_1_2_3_4_5, rel_2cl_2_2_1, FULL, rel_inner_replacement7_5_5_1, FULL, {1, 4, 5, 6, 0}));

RAM* scc4702 = new RAM(false, 7);
scc4702->add_relation(rel_4cl_4_1_2_3_4, true);
scc4702->add_relation(rel_4cl_4_4_3_2, true);
scc4702->add_rule(new parallel_acopy(rel_4cl_4_4_3_2, rel_4cl_4_1_2_3_4, DELTA, {3, 2, 1, 4, 0}));

RAM* scc4703 = new RAM(false, 11);
scc4703->add_relation(rel_graph_2_1_2, false);
scc4703->add_relation(rel_2cl_2_2_1, true);
scc4703->add_rule(new parallel_copy_filter(rel_2cl_2_2_1, rel_graph_2_1_2, FULL, {0, 1}, [](const u64* const data){ return ((s32)(data[0]&0xffffffff) < (s32)(data[1]&0xffffffff)); }));

RAM* scc4704 = new RAM(false, 15);
scc4704->add_relation(rel_2cl_2_2, true);
scc4704->add_relation(rel_2cl_2_2_1, true);
scc4704->add_rule(new parallel_acopy(rel_2cl_2_2, rel_2cl_2_2_1, DELTA, {0, 2, 1}));

RAM* scc4705 = new RAM(false, 17);
scc4705->add_relation(rel_2cl_2_1, true);
scc4705->add_relation(rel_2cl_2_2_1, true);
scc4705->add_rule(new parallel_acopy(rel_2cl_2_1, rel_2cl_2_2_1, DELTA, {1, 2, 0}));

RAM* scc4706 = new RAM(false, 2);
scc4706->add_relation(rel_inter_body26_3_1_3, false);
scc4706->add_relation(rel_3cl_3_1_2_3, true);
scc4706->add_relation(rel_2cl_2_2_1, false);
scc4706->add_rule(new parallel_join(rel_3cl_3_1_2_3, rel_2cl_2_2_1, FULL, rel_inter_body26_3_1_3, FULL, {1, 4, 0}));

RAM* scc4707 = new RAM(false, 6);
scc4707->add_relation(rel_graph_2_1_2, true);
scc4707->add_relation(rel_graph_2_2_1, true);
scc4707->add_rule(new parallel_acopy(rel_graph_2_1_2, rel_graph_2_2_1, DELTA, {1, 0, 2}));

RAM* scc4708 = new RAM(false, 10);
scc4708->add_relation(rel_inner_replacement21_4_4_1, false);
scc4708->add_relation(rel_4cl_4_1_2_3_4, true);
scc4708->add_relation(rel_2cl_2_2_1, false);
scc4708->add_rule(new parallel_join(rel_4cl_4_1_2_3_4, rel_2cl_2_2_1, FULL, rel_inner_replacement21_4_4_1, FULL, {1, 4, 5, 0}));

RAM* scc4709 = new RAM(false, 14);
scc4709->add_relation(rel_inner_replacement7_5_1_2_3_4_5, true);
scc4709->add_relation(rel_4cl_4_3_2_1, false);
scc4709->add_relation(rel_4cl_4_4_3_2, false);
scc4709->add_rule(new parallel_join(rel_inner_replacement7_5_1_2_3_4_5, rel_4cl_4_3_2_1, FULL, rel_4cl_4_4_3_2, FULL, {6, 2, 1, 0, 4}));

RAM* scc4710 = new RAM(false, 4);
scc4710->add_relation(rel_inter_body26_3_1_3, true);
scc4710->add_relation(rel_inter_body26_3_1_2_3, true);
scc4710->add_rule(new parallel_acopy(rel_inter_body26_3_1_3, rel_inter_body26_3_1_2_3, DELTA, {0, 2, 3, 1}));

RAM* scc4711 = new RAM(false, 8);
scc4711->add_relation(rel_3cl_3_3_2, false);
scc4711->add_relation(rel_3cl_3_2_1, false);
scc4711->add_relation(rel_inner_replacement21_4_1_2_3_4, true);
scc4711->add_rule(new parallel_join(rel_inner_replacement21_4_1_2_3_4, rel_3cl_3_2_1, FULL, rel_3cl_3_3_2, FULL, {5, 1, 0, 3}));

RAM* scc4712 = new RAM(false, 12);
scc4712->add_relation(rel_inner_replacement7_5_5_1, true);
scc4712->add_relation(rel_inner_replacement7_5_1_2_3_4_5, true);
scc4712->add_rule(new parallel_acopy(rel_inner_replacement7_5_5_1, rel_inner_replacement7_5_1_2_3_4_5, DELTA, {4, 0, 5, 1, 2, 3}));

RAM* scc4713 = new RAM(false, 16);
scc4713->add_relation(rel_4cl_4_1_2_3_4, true);
scc4713->add_relation(rel_4cl_4_3_2_1, true);
scc4713->add_rule(new parallel_acopy(rel_4cl_4_3_2_1, rel_4cl_4_1_2_3_4, DELTA, {2, 1, 0, 4, 3}));

LIE* lie = new LIE();
lie->add_relation(rel_4cl_4_4_3_2);
lie->add_relation(rel_2cl_2_2_1);
lie->add_relation(rel_inner_replacement21_4_1_2_3_4);
lie->add_relation(rel_4cl_4_3_2_1);
lie->add_relation(rel_4cl_4_1_2_3_4);
lie->add_relation(rel_5cl_5_1_2_3_4_5);
lie->add_relation(rel_2cl_2_1);
lie->add_relation(rel_inner_replacement21_4_4_1);
lie->add_relation(rel_2cl_2_2);
lie->add_relation(rel_graph_2_2_1);
lie->add_relation(rel_graph_2_1_2);
lie->add_relation(rel_inter_body26_3_1_2_3);
lie->add_relation(rel_3cl_3_1_2_3);
lie->add_relation(rel_inner_replacement7_5_1_2_3_4_5);
lie->add_relation(rel_3cl_3_2_1);
lie->add_relation(rel_inner_replacement7_5_5_1);
lie->add_relation(rel_inter_body26_3_1_3);
lie->add_relation(rel_3cl_3_3_2);
lie->add_scc(scc4696);
lie->add_scc(scc4697);
lie->add_scc(scc4698);
lie->add_scc(scc4699);
lie->add_scc(scc4700);
lie->add_scc(scc4701);
lie->add_scc(scc4702);
lie->add_scc(scc4703);
lie->add_scc(scc4704);
lie->add_scc(scc4705);
lie->add_scc(scc4706);
lie->add_scc(scc4707);
lie->add_scc(scc4708);
lie->add_scc(scc4709);
lie->add_scc(scc4710);
lie->add_scc(scc4711);
lie->add_scc(scc4712);
lie->add_scc(scc4713);
lie->add_scc_dependance(scc4696, scc4708);
lie->add_scc_dependance(scc4696, scc4706);
lie->add_scc_dependance(scc4696, scc4705);
lie->add_scc_dependance(scc4696, scc4704);
lie->add_scc_dependance(scc4696, scc4701);
lie->add_scc_dependance(scc4697, scc4711);
lie->add_scc_dependance(scc4698, scc4708);
lie->add_scc_dependance(scc4699, scc4711);
lie->add_scc_dependance(scc4700, scc4710);
lie->add_scc_dependance(scc4702, scc4709);
lie->add_scc_dependance(scc4703, scc4708);
lie->add_scc_dependance(scc4703, scc4706);
lie->add_scc_dependance(scc4703, scc4705);
lie->add_scc_dependance(scc4703, scc4704);
lie->add_scc_dependance(scc4703, scc4701);
lie->add_scc_dependance(scc4704, scc4700);
lie->add_scc_dependance(scc4705, scc4700);
lie->add_scc_dependance(scc4706, scc4699);
lie->add_scc_dependance(scc4706, scc4697);
lie->add_scc_dependance(scc4707, scc4703);
lie->add_scc_dependance(scc4708, scc4713);
lie->add_scc_dependance(scc4708, scc4702);
lie->add_scc_dependance(scc4709, scc4712);
lie->add_scc_dependance(scc4710, scc4706);
lie->add_scc_dependance(scc4711, scc4698);
lie->add_scc_dependance(scc4712, scc4701);
lie->add_scc_dependance(scc4713, scc4709);





    lie->set_debug_output_filename(argv[1]);
    lie->set_comm(mcomm);
    lie->set_batch_size(100);
    lie->execute();
    lie->print_all_relation_size();

    delete lie;

    mcomm.destroy();
    return 0;
}
