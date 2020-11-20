// Compilation template for slog daemon
#include "../src/parallel_RA_inc.h"

int main(int argc, char **argv)
{
  mpi_comm mcomm;
  mcomm.create(argc, argv);

relation* rel_inner_replacement53_8_6_5_4_3_2_7 = new relation(6, false, 8, 266, "rel_inner_replacement53_8_6_5_4_3_2_7", "../data/worstcase-110-terms-5-m//inner-replacement53_8_58.dat", FULL);
relation* rel_seq_2_0 = new relation(1, false, 2, 257, "rel_seq_2_0", "../data/worstcase-110-terms-5-m//seq_2_57.dat", FULL);
relation* rel_t_ret_to_ret_12_9_8_7 = new relation(3, false, 12, 262, "rel_t_ret_to_ret_12_9_8_7", "../data/worstcase-110-terms-5-m//t-ret-to-ret_12_56.dat", FULL);
relation* rel_ret_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12 = new relation(12, true, 12, 271, "rel_ret_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12", "../data/worstcase-110-terms-5-m//ret-to-ret_12_55.dat", FULL);
relation* rel_var_to_var_12_6_5_4_3_2_1 = new relation(6, false, 12, 278, "rel_var_to_var_12_6_5_4_3_2_1", "../data/worstcase-110-terms-5-m//var-to-var_12_54.dat", FULL);
relation* rel_app_2_1_2 = new relation(2, true, 2, 261, "rel_app_2_1_2", "../data/worstcase-110-terms-5-m//app_2_53.dat", FULL);
relation* rel_ret_to_ret_12_6_5_4_3_2_1 = new relation(6, false, 12, 271, "rel_ret_to_ret_12_6_5_4_3_2_1", "../data/worstcase-110-terms-5-m//ret-to-ret_12_52.dat", FULL);
relation* rel_let_3_3 = new relation(1, false, 3, 282, "rel_let_3_3", "../data/worstcase-110-terms-5-m//let_3_51.dat", FULL);
relation* rel_ref_1_0 = new relation(1, false, 1, 265, "rel_ref_1_0", "../data/worstcase-110-terms-5-m//ref_1_50.dat", FULL);
relation* rel_const_1_1 = new relation(1, true, 1, 281, "rel_const_1_1", "../data/worstcase-110-terms-5-m//const_1_49.dat", FULL);
relation* rel_app_2_1 = new relation(1, false, 2, 261, "rel_app_2_1", "../data/worstcase-110-terms-5-m//app_2_48.dat", FULL);
relation* rel_free_var_prop_11_1 = new relation(1, false, 11, 274, "rel_free_var_prop_11_1", "../data/worstcase-110-terms-5-m//free-var-prop_11_47.dat", FULL);
relation* rel_t_ret_to_var_12_12_11_10_9_8_7 = new relation(6, false, 12, 258, "rel_t_ret_to_var_12_12_11_10_9_8_7", "../data/worstcase-110-terms-5-m//t-ret-to-var_12_46.dat", FULL);
relation* rel_seq_2_1 = new relation(1, false, 2, 257, "rel_seq_2_1", "../data/worstcase-110-terms-5-m//seq_2_45.dat", FULL);
relation* rel_inter_head1420_7_1_2_3_4_5_6_7 = new relation(7, true, 7, 267, "rel_inter_head1420_7_1_2_3_4_5_6_7", "../data/worstcase-110-terms-5-m//inter-head1420_7_44.dat", FULL);
relation* rel_inter_body1436_2_1 = new relation(1, false, 2, 259, "rel_inter_body1436_2_1", "../data/worstcase-110-terms-5-m//inter-body1436_2_43.dat", FULL);
relation* rel_t_ret_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12 = new relation(12, true, 12, 262, "rel_t_ret_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12", "../data/worstcase-110-terms-5-m//t-ret-to-ret_12_42.dat", FULL);
relation* rel_lam_2_2 = new relation(1, false, 2, 270, "rel_lam_2_2", "../data/worstcase-110-terms-5-m//lam_2_41.dat", FULL);
relation* rel_seq_2_1_2 = new relation(2, true, 2, 257, "rel_seq_2_1_2", "../data/worstcase-110-terms-5-m//seq_2_40.dat", FULL);
relation* rel_program_1_1 = new relation(1, true, 1, 272, "rel_program_1_1", "../data/worstcase-110-terms-5-m//program_1_39.dat", FULL);
relation* rel_reachable_6_1_2_3_4_5_6 = new relation(6, true, 6, 280, "rel_reachable_6_1_2_3_4_5_6", "../data/worstcase-110-terms-5-m//reachable_6_38.dat", FULL);
relation* rel_ret_to_ret_9_3_2_1 = new relation(3, false, 9, 284, "rel_ret_to_ret_9_3_2_1", "../data/worstcase-110-terms-5-m//ret-to-ret_9_37.dat", FULL);
relation* rel_inner_replacement53_8_1_2_3_4_5_6_7_8 = new relation(8, true, 8, 266, "rel_inner_replacement53_8_1_2_3_4_5_6_7_8", "../data/worstcase-110-terms-5-m//inner-replacement53_8_36.dat", FULL);
relation* rel_inter_body1436_2_1_2 = new relation(2, true, 2, 259, "rel_inter_body1436_2_1_2", "../data/worstcase-110-terms-5-m//inter-body1436_2_35.dat", FULL);
relation* rel_const_1_0 = new relation(1, false, 1, 281, "rel_const_1_0", "../data/worstcase-110-terms-5-m//const_1_34.dat", FULL);
relation* rel_t_ret_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12 = new relation(12, true, 12, 258, "rel_t_ret_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12", "../data/worstcase-110-terms-5-m//t-ret-to-var_12_33.dat", FULL);
relation* rel_app_2_2 = new relation(1, false, 2, 261, "rel_app_2_2", "../data/worstcase-110-terms-5-m//app_2_32.dat", FULL);
relation* rel_let_3_2 = new relation(1, false, 3, 282, "rel_let_3_2", "../data/worstcase-110-terms-5-m//let_3_31.dat", FULL);
relation* rel_inter_head1430_15_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15 = new relation(15, true, 15, 276, "rel_inter_head1430_15_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15", "../data/worstcase-110-terms-5-m//inter-head1430_15_30.dat", FULL);
relation* rel_ret_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12 = new relation(12, true, 12, 285, "rel_ret_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12", "../data/worstcase-110-terms-5-m//ret-to-var_12_29.dat", FULL);
relation* rel_app_step_12_1_2_3_4_5_6_7_8_9_10_11_12 = new relation(12, true, 12, 260, "rel_app_step_12_1_2_3_4_5_6_7_8_9_10_11_12", "../data/worstcase-110-terms-5-m//app-step_12_28.dat", FULL);
relation* rel_ret_to_var_12_6_5_4_3_2_1 = new relation(6, false, 12, 285, "rel_ret_to_var_12_6_5_4_3_2_1", "../data/worstcase-110-terms-5-m//ret-to-var_12_27.dat", FULL);
relation* rel_var_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12 = new relation(12, true, 12, 278, "rel_var_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12", "../data/worstcase-110-terms-5-m//var-to-var_12_26.dat", FULL);
relation* rel_t_ret_to_ret_12_12_11_10_9_8_7 = new relation(6, false, 12, 262, "rel_t_ret_to_ret_12_12_11_10_9_8_7", "../data/worstcase-110-terms-5-m//t-ret-to-ret_12_25.dat", FULL);
relation* rel_inter_body1423_3_2_1 = new relation(2, false, 3, 264, "rel_inter_body1423_3_2_1", "../data/worstcase-110-terms-5-m//inter-body1423_3_24.dat", FULL);
relation* rel_lam_2_1_2 = new relation(2, true, 2, 270, "rel_lam_2_1_2", "../data/worstcase-110-terms-5-m//lam_2_23.dat", FULL);
relation* rel_inter_head1427_7_1_2_3_4_5_6_7 = new relation(7, true, 7, 273, "rel_inter_head1427_7_1_2_3_4_5_6_7", "../data/worstcase-110-terms-5-m//inter-head1427_7_22.dat", FULL);
relation* rel_inner_replacement55_13_1 = new relation(1, false, 13, 269, "rel_inner_replacement55_13_1", "../data/worstcase-110-terms-5-m//inner-replacement55_13_21.dat", FULL);
relation* rel_ret_to_ret_9_1_2_3_4_5_6_7_8_9 = new relation(9, true, 9, 284, "rel_ret_to_ret_9_1_2_3_4_5_6_7_8_9", "../data/worstcase-110-terms-5-m//ret-to-ret_9_20.dat", FULL);
relation* rel_ref_1_1 = new relation(1, true, 1, 265, "rel_ref_1_1", "../data/worstcase-110-terms-5-m//ref_1_19.dat", FULL);
relation* rel_free_var_prop_11_1_2_3_4_5_6_7_8_9_10_11 = new relation(11, true, 11, 274, "rel_free_var_prop_11_1_2_3_4_5_6_7_8_9_10_11", "../data/worstcase-110-terms-5-m//free-var-prop_11_18.dat", FULL);
relation* rel_lam_2_0 = new relation(1, false, 2, 270, "rel_lam_2_0", "../data/worstcase-110-terms-5-m//lam_2_17.dat", FULL);
relation* rel_var_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12 = new relation(12, true, 12, 256, "rel_var_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12", "../data/worstcase-110-terms-5-m//var-to-ret_12_16.dat", FULL);
relation* rel_free_2_1_2 = new relation(2, true, 2, 275, "rel_free_2_1_2", "../data/worstcase-110-terms-5-m//free_2_15.dat", FULL);
relation* rel_inner_replacement55_13_1_2_3_4_5_6_7_8_9_10_11_12_13 = new relation(13, true, 13, 269, "rel_inner_replacement55_13_1_2_3_4_5_6_7_8_9_10_11_12_13", "../data/worstcase-110-terms-5-m//inner-replacement55_13_14.dat", FULL);
relation* rel_var_to_ret_12_6_5_4_3_2_1 = new relation(6, false, 12, 256, "rel_var_to_ret_12_6_5_4_3_2_1", "../data/worstcase-110-terms-5-m//var-to-ret_12_13.dat", FULL);
relation* rel_free_2_2 = new relation(1, false, 2, 275, "rel_free_2_2", "../data/worstcase-110-terms-5-m//free_2_12.dat", FULL);
relation* rel_producer_6_6_5_4_3_2_1 = new relation(6, true, 6, 268, "rel_producer_6_6_5_4_3_2_1", "../data/worstcase-110-terms-5-m//producer_6_11.dat", FULL);
relation* rel_inter_body1416_3_3_2 = new relation(2, false, 3, 283, "rel_inter_body1416_3_3_2", "../data/worstcase-110-terms-5-m//inter-body1416_3_10.dat", FULL);
relation* rel_app_2_0 = new relation(1, false, 2, 261, "rel_app_2_0", "../data/worstcase-110-terms-5-m//app_2_9.dat", FULL);
relation* rel_inter_body1416_3_1_2_3 = new relation(3, true, 3, 283, "rel_inter_body1416_3_1_2_3", "../data/worstcase-110-terms-5-m//inter-body1416_3_8.dat", FULL);
relation* rel_reachable_6_1 = new relation(1, false, 6, 280, "rel_reachable_6_1", "../data/worstcase-110-terms-5-m//reachable_6_7.dat", FULL);
relation* rel_inner_replacement54_13_1_2_3_4_5_6_7_8_9_10_11_12_13 = new relation(13, true, 13, 277, "rel_inner_replacement54_13_1_2_3_4_5_6_7_8_9_10_11_12_13", "../data/worstcase-110-terms-5-m//inner-replacement54_13_6.dat", FULL);
relation* rel_inter_body1423_3_1_2_3 = new relation(3, true, 3, 264, "rel_inter_body1423_3_1_2_3", "../data/worstcase-110-terms-5-m//inter-body1423_3_5.dat", FULL);
relation* rel_inter_head1433_8_1_2_3_4_5_6_7_8 = new relation(8, true, 8, 263, "rel_inter_head1433_8_1_2_3_4_5_6_7_8", "../data/worstcase-110-terms-5-m//inter-head1433_8_4.dat", FULL);
relation* rel_let_3_1_2_3 = new relation(3, true, 3, 282, "rel_let_3_1_2_3", "../data/worstcase-110-terms-5-m//let_3_3.dat", FULL);
relation* rel_let_3_0 = new relation(1, false, 3, 282, "rel_let_3_0", "../data/worstcase-110-terms-5-m//let_3_2.dat", FULL);
relation* rel_inner_replacement54_13_6_5_4_3_2_1 = new relation(6, false, 13, 277, "rel_inner_replacement54_13_6_5_4_3_2_1", "../data/worstcase-110-terms-5-m//inner-replacement54_13_1.dat", FULL);
relation* rel_seq_2_2 = new relation(1, false, 2, 257, "rel_seq_2_2", "../data/worstcase-110-terms-5-m//seq_2_0.dat", FULL);

RAM* scc4991 = new RAM(false, 1);
scc4991->add_relation(rel_ref_1_1, true);
scc4991->add_relation(rel_ref_1_0, true);
scc4991->add_rule(new parallel_acopy(rel_ref_1_0, rel_ref_1_1, DELTA, {1, 0}));

RAM* scc4992 = new RAM(false, 5);
scc4992->add_relation(rel_const_1_0, true);
scc4992->add_relation(rel_const_1_1, true);
scc4992->add_rule(new parallel_acopy(rel_const_1_0, rel_const_1_1, DELTA, {1, 0}));

RAM* scc4993 = new RAM(false, 9);
scc4993->add_relation(rel_let_3_1_2_3, true);
scc4993->add_relation(rel_let_3_3, true);
scc4993->add_rule(new parallel_acopy(rel_let_3_3, rel_let_3_1_2_3, DELTA, {2, 3, 0, 1}));

RAM* scc4994 = new RAM(false, 13);
scc4994->add_relation(rel_seq_2_2, true);
scc4994->add_relation(rel_seq_2_1_2, true);
scc4994->add_rule(new parallel_acopy(rel_seq_2_2, rel_seq_2_1_2, DELTA, {1, 2, 0}));

RAM* scc4995 = new RAM(false, 18);
scc4995->add_relation(rel_let_3_0, true);
scc4995->add_relation(rel_let_3_1_2_3, true);
scc4995->add_rule(new parallel_acopy(rel_let_3_0, rel_let_3_1_2_3, DELTA, {3, 0, 1, 2}));

RAM* scc4996 = new RAM(false, 3);
scc4996->add_relation(rel_lam_2_0, true);
scc4996->add_relation(rel_lam_2_1_2, true);
scc4996->add_rule(new parallel_acopy(rel_lam_2_0, rel_lam_2_1_2, DELTA, {2, 0, 1}));

RAM* scc4997 = new RAM(false, 7);
scc4997->add_relation(rel_app_2_0, true);
scc4997->add_relation(rel_app_2_1_2, true);
scc4997->add_rule(new parallel_acopy(rel_app_2_0, rel_app_2_1_2, DELTA, {2, 0, 1}));

RAM* scc4998 = new RAM(false, 11);
scc4998->add_relation(rel_let_3_1_2_3, true);
scc4998->add_relation(rel_let_3_2, true);
scc4998->add_rule(new parallel_acopy(rel_let_3_2, rel_let_3_1_2_3, DELTA, {1, 3, 0, 2}));

RAM* scc4999 = new RAM(true, 15);
scc4999->add_relation(rel_seq_2_2, false);
scc4999->add_relation(rel_inter_body1423_3_1_2_3, true);
scc4999->add_relation(rel_inter_body1416_3_1_2_3, true);
scc4999->add_relation(rel_inter_body1416_3_3_2, true);
scc4999->add_relation(rel_free_2_2, true);
scc4999->add_relation(rel_free_2_1_2, true);
scc4999->add_relation(rel_inter_body1423_3_2_1, true);
scc4999->add_relation(rel_let_3_2, false);
scc4999->add_relation(rel_app_2_2, false);
scc4999->add_relation(rel_lam_2_2, false);
scc4999->add_relation(rel_seq_2_1, false);
scc4999->add_relation(rel_app_2_1, false);
scc4999->add_relation(rel_let_3_3, false);
scc4999->add_rule(new parallel_acopy(rel_free_2_2, rel_free_2_1_2, DELTA, {1, 2, 0}));
scc4999->add_rule(new parallel_join(rel_free_2_1_2, rel_free_2_2, DELTA, rel_app_2_2, FULL, {2, 3}));
scc4999->add_rule(new parallel_acopy(rel_inter_body1423_3_2_1, rel_inter_body1423_3_1_2_3, DELTA, {1, 0, 3, 2}));
scc4999->add_rule(new parallel_join(rel_inter_body1423_3_1_2_3, rel_lam_2_2, FULL, rel_free_2_2, DELTA, {2, 4, 1}));
scc4999->add_rule(new parallel_join(rel_free_2_1_2, rel_free_2_2, DELTA, rel_seq_2_1, FULL, {2, 3}));
scc4999->add_rule(new parallel_copy_filter(rel_free_2_1_2, rel_inter_body1416_3_3_2, DELTA, {0, 3}, [](const u64* const data){ return !(data[0] == data[1]); }));
scc4999->add_rule(new parallel_join(rel_free_2_1_2, rel_free_2_2, DELTA, rel_app_2_1, FULL, {2, 3}));
scc4999->add_rule(new parallel_acopy(rel_inter_body1416_3_3_2, rel_inter_body1416_3_1_2_3, DELTA, {2, 1, 3, 0}));
scc4999->add_rule(new parallel_join(rel_free_2_1_2, rel_seq_2_2, FULL, rel_free_2_2, DELTA, {4, 1}));
scc4999->add_rule(new parallel_copy_filter(rel_free_2_1_2, rel_inter_body1423_3_2_1, DELTA, {0, 3}, [](const u64* const data){ return !(data[0] == data[1]); }));
scc4999->add_rule(new parallel_join(rel_free_2_1_2, rel_free_2_2, DELTA, rel_let_3_2, FULL, {2, 3}));
scc4999->add_rule(new parallel_join(rel_inter_body1416_3_1_2_3, rel_free_2_2, DELTA, rel_let_3_3, FULL, {3, 4, 2}));

RAM* scc5000 = new RAM(false, 20);
scc5000->add_relation(rel_ret_to_ret_9_1_2_3_4_5_6_7_8_9, true);
scc5000->add_relation(rel_ret_to_ret_9_3_2_1, true);
scc5000->add_rule(new parallel_acopy(rel_ret_to_ret_9_3_2_1, rel_ret_to_ret_9_1_2_3_4_5_6_7_8_9, DELTA, {2, 1, 0, 9, 3, 4, 5, 6, 7, 8}));

RAM* scc5001 = new RAM(false, 17);
scc5001->add_relation(rel_free_2_2, false);
scc5001->add_relation(rel_lam_2_0, false);
scc5001->add_relation(rel_inter_body1436_2_1_2, true);
scc5001->add_rule(new parallel_join(rel_inter_body1436_2_1_2, rel_lam_2_0, FULL, rel_free_2_2, FULL, {0, 4}));

RAM* scc5002 = new RAM(false, 21);
scc5002->add_relation(rel_app_2_1, true);
scc5002->add_relation(rel_app_2_1_2, true);
scc5002->add_rule(new parallel_acopy(rel_app_2_1, rel_app_2_1_2, DELTA, {0, 2, 1}));

RAM* scc5003 = new RAM(false, 2);
scc5003->add_relation(rel_app_step_12_1_2_3_4_5_6_7_8_9_10_11_12, true);
scc5003->add_relation(rel_inter_head1430_15_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15, false);
scc5003->add_rule(new parallel_copy(rel_app_step_12_1_2_3_4_5_6_7_8_9_10_11_12, rel_inter_head1430_15_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15, FULL, {12, 13, 9, 2, 10, 11, 14, 12, 13, 9, 2, 10}));

RAM* scc5004 = new RAM(false, 6);
scc5004->add_relation(rel_app_2_2, true);
scc5004->add_relation(rel_app_2_1_2, true);
scc5004->add_rule(new parallel_acopy(rel_app_2_2, rel_app_2_1_2, DELTA, {1, 2, 0}));

RAM* scc5005 = new RAM(false, 10);
scc5005->add_relation(rel_seq_2_1_2, true);
scc5005->add_relation(rel_seq_2_0, true);
scc5005->add_rule(new parallel_acopy(rel_seq_2_0, rel_seq_2_1_2, DELTA, {2, 0, 1}));

RAM* scc5006 = new RAM(false, 14);
scc5006->add_relation(rel_lam_2_1_2, true);
scc5006->add_relation(rel_lam_2_2, true);
scc5006->add_rule(new parallel_acopy(rel_lam_2_2, rel_lam_2_1_2, DELTA, {1, 2, 0}));

RAM* scc5007 = new RAM(false, 19);
scc5007->add_relation(rel_inter_body1436_2_1_2, true);
scc5007->add_relation(rel_inter_body1436_2_1, true);
scc5007->add_rule(new parallel_acopy(rel_inter_body1436_2_1, rel_inter_body1436_2_1_2, DELTA, {0, 2, 1}));

RAM* scc5008 = new RAM(false, 4);
scc5008->add_relation(rel_reachable_6_1_2_3_4_5_6, true);
scc5008->add_relation(rel_program_1_1, false);
scc5008->add_rule(new parallel_copy(rel_reachable_6_1_2_3_4_5_6, rel_program_1_1, FULL, {0, 0, 0, 0, 0, 0}));

RAM* scc5009 = new RAM(false, 8);
scc5009->add_relation(rel_free_2_1_2, true);
scc5009->add_relation(rel_ref_1_1, false);
scc5009->add_rule(new parallel_copy(rel_free_2_1_2, rel_ref_1_1, FULL, {0, 1}));

RAM* scc5010 = new RAM(true, 12);
scc5010->add_relation(rel_inner_replacement54_13_6_5_4_3_2_1, true);
scc5010->add_relation(rel_let_3_0, false);
scc5010->add_relation(rel_inter_head1433_8_1_2_3_4_5_6_7_8, true);
scc5010->add_relation(rel_inner_replacement54_13_1_2_3_4_5_6_7_8_9_10_11_12_13, true);
scc5010->add_relation(rel_reachable_6_1, true);
scc5010->add_relation(rel_app_2_0, false);
scc5010->add_relation(rel_producer_6_6_5_4_3_2_1, true);
scc5010->add_relation(rel_var_to_ret_12_6_5_4_3_2_1, true);
scc5010->add_relation(rel_inner_replacement55_13_1_2_3_4_5_6_7_8_9_10_11_12_13, true);
scc5010->add_relation(rel_var_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12, true);
scc5010->add_relation(rel_lam_2_0, false);
scc5010->add_relation(rel_free_var_prop_11_1_2_3_4_5_6_7_8_9_10_11, true);
scc5010->add_relation(rel_inner_replacement55_13_1, true);
scc5010->add_relation(rel_inter_head1427_7_1_2_3_4_5_6_7, true);
scc5010->add_relation(rel_t_ret_to_ret_12_12_11_10_9_8_7, true);
scc5010->add_relation(rel_var_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12, true);
scc5010->add_relation(rel_ret_to_var_12_6_5_4_3_2_1, true);
scc5010->add_relation(rel_ret_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12, true);
scc5010->add_relation(rel_inter_head1430_15_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15, true);
scc5010->add_relation(rel_t_ret_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12, true);
scc5010->add_relation(rel_const_1_0, false);
scc5010->add_relation(rel_inner_replacement53_8_1_2_3_4_5_6_7_8, true);
scc5010->add_relation(rel_ret_to_ret_9_3_2_1, false);
scc5010->add_relation(rel_reachable_6_1_2_3_4_5_6, true);
scc5010->add_relation(rel_t_ret_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12, true);
scc5010->add_relation(rel_inter_body1436_2_1, false);
scc5010->add_relation(rel_inter_head1420_7_1_2_3_4_5_6_7, true);
scc5010->add_relation(rel_t_ret_to_var_12_12_11_10_9_8_7, true);
scc5010->add_relation(rel_free_var_prop_11_1, true);
scc5010->add_relation(rel_ref_1_0, false);
scc5010->add_relation(rel_ret_to_ret_12_6_5_4_3_2_1, true);
scc5010->add_relation(rel_var_to_var_12_6_5_4_3_2_1, true);
scc5010->add_relation(rel_ret_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12, true);
scc5010->add_relation(rel_t_ret_to_ret_12_9_8_7, true);
scc5010->add_relation(rel_seq_2_0, false);
scc5010->add_relation(rel_inner_replacement53_8_6_5_4_3_2_7, true);
scc5010->add_rule(new parallel_join(rel_t_ret_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12, rel_var_to_ret_12_6_5_4_3_2_1, FULL, rel_t_ret_to_var_12_12_11_10_9_8_7, DELTA, {14, 15, 16, 17, 18, 19, 7, 8, 9, 10, 11, 12}));
scc5010->add_rule(new parallel_join(rel_inter_head1433_8_1_2_3_4_5_6_7_8, rel_let_3_0, FULL, rel_reachable_6_1, DELTA, {7, 2, 1, 6, 8, 9, 5, 3}));
scc5010->add_rule(new parallel_acopy(rel_t_ret_to_ret_12_9_8_7, rel_t_ret_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12, DELTA, {8, 7, 6, 12, 0, 1, 2, 3, 4, 5, 9, 10, 11}));
scc5010->add_rule(new parallel_join(rel_inner_replacement55_13_1_2_3_4_5_6_7_8_9_10_11_12_13, rel_producer_6_6_5_4_3_2_1, DELTA, rel_inner_replacement54_13_6_5_4_3_2_1, DELTA, {5, 4, 3, 2, 1, 0, 8, 9, 10, 11, 12, 13, 14}));
scc5010->add_rule(new parallel_acopy(rel_var_to_ret_12_6_5_4_3_2_1, rel_var_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12, DELTA, {5, 4, 3, 2, 1, 0, 12, 6, 7, 8, 9, 10, 11}));
scc5010->add_rule(new parallel_acopy(rel_free_var_prop_11_1, rel_free_var_prop_11_1_2_3_4_5_6_7_8_9_10_11, DELTA, {0, 11, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}));
scc5010->add_rule(new parallel_acopy(rel_inner_replacement53_8_6_5_4_3_2_7, rel_inner_replacement53_8_1_2_3_4_5_6_7_8, DELTA, {5, 4, 3, 2, 1, 6, 8, 0, 7}));
scc5010->add_rule(new parallel_acopy(rel_var_to_var_12_6_5_4_3_2_1, rel_var_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12, DELTA, {5, 4, 3, 2, 1, 0, 12, 6, 7, 8, 9, 10, 11}));
scc5010->add_rule(new parallel_copy(rel_reachable_6_1_2_3_4_5_6, rel_inter_head1420_7_1_2_3_4_5_6_7, DELTA, {1, 6, 3, 2, 4, 5}));
scc5010->add_rule(new parallel_copy(rel_reachable_6_1_2_3_4_5_6, rel_inter_head1433_8_1_2_3_4_5_6_7_8, DELTA, {1, 6, 3, 0, 4, 5}));
scc5010->add_rule(new parallel_acopy(rel_t_ret_to_var_12_12_11_10_9_8_7, rel_t_ret_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12, DELTA, {11, 10, 9, 8, 7, 6, 12, 0, 1, 2, 3, 4, 5}));
scc5010->add_rule(new parallel_join(rel_t_ret_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12, rel_producer_6_6_5_4_3_2_1, FULL, rel_ret_to_ret_12_6_5_4_3_2_1, DELTA, {5, 4, 3, 2, 1, 0, 8, 9, 10, 11, 12, 13}));
scc5010->add_rule(new parallel_join(rel_t_ret_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12, rel_t_ret_to_ret_12_12_11_10_9_8_7, FULL, rel_ret_to_var_12_6_5_4_3_2_1, DELTA, {7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19}));
scc5010->add_rule(new parallel_copy(rel_ret_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12, rel_inter_head1430_15_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15, DELTA, {0, 13, 9, 2, 10, 11, 8, 12, 13, 9, 2, 10}));
scc5010->add_rule(new parallel_acopy(rel_inner_replacement54_13_6_5_4_3_2_1, rel_inner_replacement54_13_1_2_3_4_5_6_7_8_9_10_11_12_13, DELTA, {5, 4, 3, 2, 1, 0, 13, 6, 7, 8, 9, 10, 11, 12}));
scc5010->add_rule(new parallel_join(rel_t_ret_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12, rel_producer_6_6_5_4_3_2_1, FULL, rel_ret_to_var_12_6_5_4_3_2_1, DELTA, {5, 4, 3, 2, 1, 0, 8, 9, 10, 11, 12, 13}));
scc5010->add_rule(new parallel_join(rel_inner_replacement53_8_1_2_3_4_5_6_7_8, rel_app_2_0, FULL, rel_reachable_6_1, DELTA, {0, 4, 5, 6, 7, 8, 1, 2}));
scc5010->add_rule(new parallel_join(rel_inner_replacement54_13_1_2_3_4_5_6_7_8_9_10_11_12_13, rel_inner_replacement53_8_6_5_4_3_2_7, DELTA, rel_t_ret_to_ret_12_12_11_10_9_8_7, FULL, {10, 11, 12, 13, 14, 15, 4, 3, 2, 1, 0, 7, 8}));
scc5010->add_rule(new parallel_acopy(rel_ret_to_var_12_6_5_4_3_2_1, rel_ret_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12, DELTA, {5, 4, 3, 2, 1, 0, 12, 6, 7, 8, 9, 10, 11}));
scc5010->add_rule(new parallel_join(rel_producer_6_6_5_4_3_2_1, rel_const_1_0, FULL, rel_reachable_6_1, DELTA, {7, 6, 5, 4, 3, 0}));
scc5010->add_rule(new parallel_acopy(rel_reachable_6_1, rel_reachable_6_1_2_3_4_5_6, DELTA, {0, 6, 1, 2, 3, 4, 5}));
scc5010->add_rule(new parallel_join(rel_inter_head1420_7_1_2_3_4_5_6_7, rel_seq_2_0, FULL, rel_reachable_6_1, DELTA, {1, 2, 6, 5, 7, 8, 4}));
scc5010->add_rule(new parallel_acopy(rel_t_ret_to_ret_12_12_11_10_9_8_7, rel_t_ret_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12, DELTA, {11, 10, 9, 8, 7, 6, 12, 0, 1, 2, 3, 4, 5}));
scc5010->add_rule(new parallel_acopy(rel_ret_to_ret_12_6_5_4_3_2_1, rel_ret_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12, DELTA, {5, 4, 3, 2, 1, 0, 12, 6, 7, 8, 9, 10, 11}));
scc5010->add_rule(new parallel_copy(rel_ret_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12, rel_producer_6_6_5_4_3_2_1, DELTA, {5, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0}));
scc5010->add_rule(new parallel_join(rel_t_ret_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12, rel_t_ret_to_var_12_12_11_10_9_8_7, FULL, rel_var_to_var_12_6_5_4_3_2_1, DELTA, {7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19}));
scc5010->add_rule(new parallel_copy(rel_free_var_prop_11_1_2_3_4_5_6_7_8_9_10_11, rel_inter_head1430_15_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15, DELTA, {1, 5, 3, 7, 6, 4, 12, 13, 9, 2, 10}));
scc5010->add_rule(new parallel_copy(rel_ret_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12, rel_inter_head1430_15_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15, DELTA, {14, 12, 13, 9, 2, 10, 12, 13, 9, 2, 10, 11}));
scc5010->add_rule(new parallel_join(rel_t_ret_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12, rel_producer_6_6_5_4_3_2_1, DELTA, rel_ret_to_ret_12_6_5_4_3_2_1, DELTA, {5, 4, 3, 2, 1, 0, 8, 9, 10, 11, 12, 13}));
scc5010->add_rule(new parallel_join(rel_t_ret_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12, rel_producer_6_6_5_4_3_2_1, DELTA, rel_ret_to_var_12_6_5_4_3_2_1, DELTA, {5, 4, 3, 2, 1, 0, 8, 9, 10, 11, 12, 13}));
scc5010->add_rule(new parallel_copy(rel_reachable_6_1_2_3_4_5_6, rel_inter_head1430_15_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15, DELTA, {14, 12, 13, 9, 2, 10}));
scc5010->add_rule(new parallel_join(rel_var_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12, rel_ref_1_0, FULL, rel_reachable_6_1, DELTA, {1, 3, 4, 5, 6, 7, 0, 3, 4, 5, 6, 7}));
scc5010->add_rule(new parallel_join(rel_producer_6_6_5_4_3_2_1, rel_lam_2_0, FULL, rel_reachable_6_1, DELTA, {8, 7, 6, 5, 4, 0}));
scc5010->add_rule(new parallel_join(rel_inner_replacement54_13_1_2_3_4_5_6_7_8_9_10_11_12_13, rel_inner_replacement53_8_6_5_4_3_2_7, FULL, rel_t_ret_to_ret_12_12_11_10_9_8_7, DELTA, {10, 11, 12, 13, 14, 15, 4, 3, 2, 1, 0, 7, 8}));
scc5010->add_rule(new parallel_copy(rel_reachable_6_1_2_3_4_5_6, rel_inter_head1427_7_1_2_3_4_5_6_7, DELTA, {0, 6, 3, 1, 4, 5}));
scc5010->add_rule(new parallel_join(rel_inner_replacement54_13_1_2_3_4_5_6_7_8_9_10_11_12_13, rel_inner_replacement53_8_6_5_4_3_2_7, DELTA, rel_t_ret_to_ret_12_12_11_10_9_8_7, DELTA, {10, 11, 12, 13, 14, 15, 4, 3, 2, 1, 0, 7, 8}));
scc5010->add_rule(new parallel_copy(rel_ret_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12, rel_inter_head1433_8_1_2_3_4_5_6_7_8, DELTA, {1, 6, 3, 0, 4, 5, 2, 6, 3, 0, 4, 5}));
scc5010->add_rule(new parallel_join(rel_var_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12, rel_inter_body1436_2_1, FULL, rel_free_var_prop_11_1, DELTA, {2, 4, 5, 6, 7, 8, 2, 9, 10, 11, 12, 13}));
scc5010->add_rule(new parallel_join(rel_t_ret_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12, rel_ret_to_ret_9_3_2_1, FULL, rel_t_ret_to_ret_12_9_8_7, DELTA, {11, 12, 13, 14, 15, 16, 4, 5, 6, 7, 8, 9}));
scc5010->add_rule(new parallel_copy(rel_reachable_6_1_2_3_4_5_6, rel_inter_head1420_7_1_2_3_4_5_6_7, DELTA, {0, 6, 3, 2, 4, 5}));
scc5010->add_rule(new parallel_join(rel_t_ret_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12, rel_producer_6_6_5_4_3_2_1, DELTA, rel_ret_to_var_12_6_5_4_3_2_1, FULL, {5, 4, 3, 2, 1, 0, 8, 9, 10, 11, 12, 13}));
scc5010->add_rule(new parallel_copy(rel_reachable_6_1_2_3_4_5_6, rel_inter_head1433_8_1_2_3_4_5_6_7_8, DELTA, {7, 6, 3, 0, 4, 5}));
scc5010->add_rule(new parallel_join(rel_t_ret_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12, rel_t_ret_to_var_12_12_11_10_9_8_7, DELTA, rel_var_to_var_12_6_5_4_3_2_1, FULL, {7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19}));
scc5010->add_rule(new parallel_join(rel_t_ret_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12, rel_t_ret_to_var_12_12_11_10_9_8_7, DELTA, rel_var_to_var_12_6_5_4_3_2_1, DELTA, {7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19}));
scc5010->add_rule(new parallel_join(rel_inner_replacement55_13_1_2_3_4_5_6_7_8_9_10_11_12_13, rel_producer_6_6_5_4_3_2_1, FULL, rel_inner_replacement54_13_6_5_4_3_2_1, DELTA, {5, 4, 3, 2, 1, 0, 8, 9, 10, 11, 12, 13, 14}));
scc5010->add_rule(new parallel_join(rel_inner_replacement55_13_1_2_3_4_5_6_7_8_9_10_11_12_13, rel_producer_6_6_5_4_3_2_1, DELTA, rel_inner_replacement54_13_6_5_4_3_2_1, FULL, {5, 4, 3, 2, 1, 0, 8, 9, 10, 11, 12, 13, 14}));
scc5010->add_rule(new parallel_join(rel_t_ret_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12, rel_var_to_ret_12_6_5_4_3_2_1, DELTA, rel_t_ret_to_var_12_12_11_10_9_8_7, FULL, {14, 15, 16, 17, 18, 19, 7, 8, 9, 10, 11, 12}));
scc5010->add_rule(new parallel_join(rel_t_ret_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12, rel_var_to_ret_12_6_5_4_3_2_1, DELTA, rel_t_ret_to_var_12_12_11_10_9_8_7, DELTA, {14, 15, 16, 17, 18, 19, 7, 8, 9, 10, 11, 12}));
scc5010->add_rule(new parallel_acopy(rel_inner_replacement55_13_1, rel_inner_replacement55_13_1_2_3_4_5_6_7_8_9_10_11_12_13, DELTA, {0, 13, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}));
scc5010->add_rule(new parallel_join(rel_t_ret_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12, rel_producer_6_6_5_4_3_2_1, DELTA, rel_ret_to_ret_12_6_5_4_3_2_1, FULL, {5, 4, 3, 2, 1, 0, 8, 9, 10, 11, 12, 13}));
scc5010->add_rule(new parallel_copy(rel_reachable_6_1_2_3_4_5_6, rel_inter_head1427_7_1_2_3_4_5_6_7, DELTA, {2, 6, 3, 1, 4, 5}));
scc5010->add_rule(new parallel_join(rel_inter_head1427_7_1_2_3_4_5_6_7, rel_app_2_0, FULL, rel_reachable_6_1, DELTA, {2, 6, 1, 5, 7, 8, 4}));
scc5010->add_rule(new parallel_join(rel_inter_head1430_15_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15, rel_lam_2_0, FULL, rel_inner_replacement55_13_1, DELTA, {15, 0, 11, 5, 8, 4, 7, 6, 1, 10, 12, 13, 14, 9, 2}));
scc5010->add_rule(new parallel_join(rel_t_ret_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12, rel_t_ret_to_ret_12_12_11_10_9_8_7, DELTA, rel_ret_to_var_12_6_5_4_3_2_1, FULL, {7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19}));
scc5010->add_rule(new parallel_join(rel_t_ret_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12, rel_t_ret_to_ret_12_12_11_10_9_8_7, DELTA, rel_ret_to_var_12_6_5_4_3_2_1, DELTA, {7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19}));

RAM* scc5011 = new RAM(false, 16);
scc5011->add_relation(rel_seq_2_1_2, true);
scc5011->add_relation(rel_seq_2_1, true);
scc5011->add_rule(new parallel_acopy(rel_seq_2_1, rel_seq_2_1_2, DELTA, {0, 2, 1}));

LIE* lie = new LIE();
lie->add_relation(rel_inner_replacement53_8_6_5_4_3_2_7);
lie->add_relation(rel_seq_2_0);
lie->add_relation(rel_t_ret_to_ret_12_9_8_7);
lie->add_relation(rel_ret_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12);
lie->add_relation(rel_var_to_var_12_6_5_4_3_2_1);
lie->add_relation(rel_app_2_1_2);
lie->add_relation(rel_ret_to_ret_12_6_5_4_3_2_1);
lie->add_relation(rel_let_3_3);
lie->add_relation(rel_ref_1_0);
lie->add_relation(rel_const_1_1);
lie->add_relation(rel_app_2_1);
lie->add_relation(rel_free_var_prop_11_1);
lie->add_relation(rel_t_ret_to_var_12_12_11_10_9_8_7);
lie->add_relation(rel_seq_2_1);
lie->add_relation(rel_inter_head1420_7_1_2_3_4_5_6_7);
lie->add_relation(rel_inter_body1436_2_1);
lie->add_relation(rel_t_ret_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12);
lie->add_relation(rel_lam_2_2);
lie->add_relation(rel_seq_2_1_2);
lie->add_relation(rel_program_1_1);
lie->add_relation(rel_reachable_6_1_2_3_4_5_6);
lie->add_relation(rel_ret_to_ret_9_3_2_1);
lie->add_relation(rel_inner_replacement53_8_1_2_3_4_5_6_7_8);
lie->add_relation(rel_inter_body1436_2_1_2);
lie->add_relation(rel_const_1_0);
lie->add_relation(rel_t_ret_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12);
lie->add_relation(rel_app_2_2);
lie->add_relation(rel_let_3_2);
lie->add_relation(rel_inter_head1430_15_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15);
lie->add_relation(rel_ret_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12);
lie->add_relation(rel_app_step_12_1_2_3_4_5_6_7_8_9_10_11_12);
lie->add_relation(rel_ret_to_var_12_6_5_4_3_2_1);
lie->add_relation(rel_var_to_var_12_1_2_3_4_5_6_7_8_9_10_11_12);
lie->add_relation(rel_t_ret_to_ret_12_12_11_10_9_8_7);
lie->add_relation(rel_inter_body1423_3_2_1);
lie->add_relation(rel_lam_2_1_2);
lie->add_relation(rel_inter_head1427_7_1_2_3_4_5_6_7);
lie->add_relation(rel_inner_replacement55_13_1);
lie->add_relation(rel_ret_to_ret_9_1_2_3_4_5_6_7_8_9);
lie->add_relation(rel_ref_1_1);
lie->add_relation(rel_free_var_prop_11_1_2_3_4_5_6_7_8_9_10_11);
lie->add_relation(rel_lam_2_0);
lie->add_relation(rel_var_to_ret_12_1_2_3_4_5_6_7_8_9_10_11_12);
lie->add_relation(rel_free_2_1_2);
lie->add_relation(rel_inner_replacement55_13_1_2_3_4_5_6_7_8_9_10_11_12_13);
lie->add_relation(rel_var_to_ret_12_6_5_4_3_2_1);
lie->add_relation(rel_free_2_2);
lie->add_relation(rel_producer_6_6_5_4_3_2_1);
lie->add_relation(rel_inter_body1416_3_3_2);
lie->add_relation(rel_app_2_0);
lie->add_relation(rel_inter_body1416_3_1_2_3);
lie->add_relation(rel_reachable_6_1);
lie->add_relation(rel_inner_replacement54_13_1_2_3_4_5_6_7_8_9_10_11_12_13);
lie->add_relation(rel_inter_body1423_3_1_2_3);
lie->add_relation(rel_inter_head1433_8_1_2_3_4_5_6_7_8);
lie->add_relation(rel_let_3_1_2_3);
lie->add_relation(rel_let_3_0);
lie->add_relation(rel_inner_replacement54_13_6_5_4_3_2_1);
lie->add_relation(rel_seq_2_2);
lie->add_scc(scc4991);
lie->add_scc(scc4992);
lie->add_scc(scc4993);
lie->add_scc(scc4994);
lie->add_scc(scc4995);
lie->add_scc(scc4996);
lie->add_scc(scc4997);
lie->add_scc(scc4998);
lie->add_scc(scc4999);
lie->add_scc(scc5000);
lie->add_scc(scc5001);
lie->add_scc(scc5002);
lie->add_scc(scc5003);
lie->add_scc(scc5004);
lie->add_scc(scc5005);
lie->add_scc(scc5006);
lie->add_scc(scc5007);
lie->add_scc(scc5008);
lie->add_scc(scc5009);
lie->add_scc(scc5010);
lie->add_scc(scc5011);
lie->add_scc_dependance(scc4991, scc5010);
lie->add_scc_dependance(scc4992, scc5010);
lie->add_scc_dependance(scc4993, scc4999);
lie->add_scc_dependance(scc4994, scc4999);
lie->add_scc_dependance(scc4995, scc5010);
lie->add_scc_dependance(scc4996, scc5010);
lie->add_scc_dependance(scc4996, scc5001);
lie->add_scc_dependance(scc4997, scc5010);
lie->add_scc_dependance(scc4998, scc4999);
lie->add_scc_dependance(scc4999, scc5001);
lie->add_scc_dependance(scc5000, scc5010);
lie->add_scc_dependance(scc5001, scc5007);
lie->add_scc_dependance(scc5002, scc4999);
lie->add_scc_dependance(scc5004, scc4999);
lie->add_scc_dependance(scc5005, scc5010);
lie->add_scc_dependance(scc5006, scc4999);
lie->add_scc_dependance(scc5007, scc5010);
lie->add_scc_dependance(scc5008, scc5010);
lie->add_scc_dependance(scc5009, scc4999);
lie->add_scc_dependance(scc5010, scc5003);
lie->add_scc_dependance(scc5011, scc4999);



  // Enable IO
 //lie->enable_share_io();
 //lie->set_output_dir(../data/worstcase-110-terms-5-m-output); // Write to this directory
 //lie->enable_share_io();
 lie->set_comm(mcomm);
 lie->set_batch_size(10);
 lie->execute();
 lie->print_all_relation_size(); // Continuously print relation sizes

 /*
 // Write to disc
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
   */
 
 delete lie;
 mcomm.destroy();
 return 0;
}
