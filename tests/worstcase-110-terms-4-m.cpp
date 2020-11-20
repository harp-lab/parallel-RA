// Compilation template for slog daemon
#include "../src/parallel_RA_inc.h"

int main(int argc, char **argv)
{
  mpi_comm mcomm;
  mcomm.create(argc, argv);

relation* rel_seq_2_0 = new relation(1, false, 2, 258, "rel_seq_2_0", "../data/worstcase-110-terms-4-m//seq_2_59.dat", FULL);
relation* rel_ret_to_var_10_1_2_3_4_5_6_7_8_9_10 = new relation(10, true, 10, 270, "rel_ret_to_var_10_1_2_3_4_5_6_7_8_9_10", "../data/worstcase-110-terms-4-m//ret-to-var_10_58.dat", FULL);
relation* rel_app_2_1_2 = new relation(2, true, 2, 259, "rel_app_2_1_2", "../data/worstcase-110-terms-4-m//app_2_57.dat", FULL);
relation* rel_inter_body1424_3_2_3 = new relation(2, false, 3, 281, "rel_inter_body1424_3_2_3", "../data/worstcase-110-terms-4-m//inter-body1424_3_56.dat", FULL);
relation* rel_let_3_3 = new relation(1, false, 3, 282, "rel_let_3_3", "../data/worstcase-110-terms-4-m//let_3_55.dat", FULL);
relation* rel_ref_1_0 = new relation(1, false, 1, 261, "rel_ref_1_0", "../data/worstcase-110-terms-4-m//ref_1_54.dat", FULL);
relation* rel_inter_head1431_6_1_2_3_4_5_6 = new relation(6, true, 6, 257, "rel_inter_head1431_6_1_2_3_4_5_6", "../data/worstcase-110-terms-4-m//inter-head1431_6_53.dat", FULL);
relation* rel_reachable_5_1 = new relation(1, false, 5, 280, "rel_reachable_5_1", "../data/worstcase-110-terms-4-m//reachable_5_52.dat", FULL);
relation* rel_var_to_ret_10_5_4_3_2_1 = new relation(5, false, 10, 283, "rel_var_to_ret_10_5_4_3_2_1", "../data/worstcase-110-terms-4-m//var-to-ret_10_51.dat", FULL);
relation* rel_const_1_1 = new relation(1, true, 1, 278, "rel_const_1_1", "../data/worstcase-110-terms-4-m//const_1_50.dat", FULL);
relation* rel_inner_replacement53_7_1_2_3_4_5_6_7 = new relation(7, true, 7, 276, "rel_inner_replacement53_7_1_2_3_4_5_6_7", "../data/worstcase-110-terms-4-m//inner-replacement53_7_49.dat", FULL);
relation* rel_inner_replacement55_11_1 = new relation(1, false, 11, 262, "rel_inner_replacement55_11_1", "../data/worstcase-110-terms-4-m//inner-replacement55_11_48.dat", FULL);
relation* rel_var_to_var_10_1_2_3_4_5_6_7_8_9_10 = new relation(10, true, 10, 274, "rel_var_to_var_10_1_2_3_4_5_6_7_8_9_10", "../data/worstcase-110-terms-4-m//var-to-var_10_47.dat", FULL);
relation* rel_app_2_1 = new relation(1, false, 2, 259, "rel_app_2_1", "../data/worstcase-110-terms-4-m//app_2_46.dat", FULL);
relation* rel_free_var_prop_9_1 = new relation(1, false, 9, 266, "rel_free_var_prop_9_1", "../data/worstcase-110-terms-4-m//free-var-prop_9_45.dat", FULL);
relation* rel_inter_body1424_3_1_2_3 = new relation(3, true, 3, 281, "rel_inter_body1424_3_1_2_3", "../data/worstcase-110-terms-4-m//inter-body1424_3_44.dat", FULL);
relation* rel_seq_2_1 = new relation(1, false, 2, 258, "rel_seq_2_1", "../data/worstcase-110-terms-4-m//seq_2_43.dat", FULL);
relation* rel_lam_2_2 = new relation(1, false, 2, 263, "rel_lam_2_2", "../data/worstcase-110-terms-4-m//lam_2_42.dat", FULL);
relation* rel_seq_2_1_2 = new relation(2, true, 2, 258, "rel_seq_2_1_2", "../data/worstcase-110-terms-4-m//seq_2_41.dat", FULL);
relation* rel_program_1_1 = new relation(1, true, 1, 264, "rel_program_1_1", "../data/worstcase-110-terms-4-m//program_1_40.dat", FULL);
relation* rel_inner_replacement54_11_5_4_3_2_1 = new relation(5, false, 11, 273, "rel_inner_replacement54_11_5_4_3_2_1", "../data/worstcase-110-terms-4-m//inner-replacement54_11_39.dat", FULL);
relation* rel_ret_to_ret_10_5_4_3_2_1 = new relation(5, false, 10, 267, "rel_ret_to_ret_10_5_4_3_2_1", "../data/worstcase-110-terms-4-m//ret-to-ret_10_38.dat", FULL);
relation* rel_var_to_var_10_5_4_3_2_1 = new relation(5, false, 10, 274, "rel_var_to_var_10_5_4_3_2_1", "../data/worstcase-110-terms-4-m//var-to-var_10_37.dat", FULL);
relation* rel_const_1_0 = new relation(1, false, 1, 278, "rel_const_1_0", "../data/worstcase-110-terms-4-m//const_1_36.dat", FULL);
relation* rel_ret_to_ret_8_3_2_1 = new relation(3, false, 8, 269, "rel_ret_to_ret_8_3_2_1", "../data/worstcase-110-terms-4-m//ret-to-ret_8_35.dat", FULL);
relation* rel_ret_to_var_10_5_4_3_2_1 = new relation(5, false, 10, 270, "rel_ret_to_var_10_5_4_3_2_1", "../data/worstcase-110-terms-4-m//ret-to-var_10_34.dat", FULL);
relation* rel_producer_5_5_4_3_2_1 = new relation(5, true, 5, 279, "rel_producer_5_5_4_3_2_1", "../data/worstcase-110-terms-4-m//producer_5_33.dat", FULL);
relation* rel_app_2_2 = new relation(1, false, 2, 259, "rel_app_2_2", "../data/worstcase-110-terms-4-m//app_2_32.dat", FULL);
relation* rel_inner_replacement53_7_6_5_4_3_2 = new relation(5, false, 7, 276, "rel_inner_replacement53_7_6_5_4_3_2", "../data/worstcase-110-terms-4-m//inner-replacement53_7_31.dat", FULL);
relation* rel_let_3_2 = new relation(1, false, 3, 282, "rel_let_3_2", "../data/worstcase-110-terms-4-m//let_3_30.dat", FULL);
relation* rel_var_to_ret_10_1_2_3_4_5_6_7_8_9_10 = new relation(10, true, 10, 283, "rel_var_to_ret_10_1_2_3_4_5_6_7_8_9_10", "../data/worstcase-110-terms-4-m//var-to-ret_10_29.dat", FULL);
relation* rel_ret_to_ret_10_1_2_3_4_5_6_7_8_9_10 = new relation(10, true, 10, 267, "rel_ret_to_ret_10_1_2_3_4_5_6_7_8_9_10", "../data/worstcase-110-terms-4-m//ret-to-ret_10_28.dat", FULL);
relation* rel_ret_to_ret_8_1_2_3_4_5_6_7_8 = new relation(8, true, 8, 269, "rel_ret_to_ret_8_1_2_3_4_5_6_7_8", "../data/worstcase-110-terms-4-m//ret-to-ret_8_27.dat", FULL);
relation* rel_free_var_prop_9_1_2_3_4_5_6_7_8_9 = new relation(9, true, 9, 266, "rel_free_var_prop_9_1_2_3_4_5_6_7_8_9", "../data/worstcase-110-terms-4-m//free-var-prop_9_26.dat", FULL);
relation* rel_lam_2_1_2 = new relation(2, true, 2, 263, "rel_lam_2_1_2", "../data/worstcase-110-terms-4-m//lam_2_25.dat", FULL);
relation* rel_inner_replacement55_11_1_2_3_4_5_6_7_8_9_10_11 = new relation(11, true, 11, 262, "rel_inner_replacement55_11_1_2_3_4_5_6_7_8_9_10_11", "../data/worstcase-110-terms-4-m//inner-replacement55_11_24.dat", FULL);
relation* rel_app_step_10_1_2_3_4_5_6_7_8_9_10 = new relation(10, true, 10, 256, "rel_app_step_10_1_2_3_4_5_6_7_8_9_10", "../data/worstcase-110-terms-4-m//app-step_10_23.dat", FULL);
relation* rel_ref_1_1 = new relation(1, true, 1, 261, "rel_ref_1_1", "../data/worstcase-110-terms-4-m//ref_1_22.dat", FULL);
relation* rel_inner_replacement54_11_1_2_3_4_5_6_7_8_9_10_11 = new relation(11, true, 11, 273, "rel_inner_replacement54_11_1_2_3_4_5_6_7_8_9_10_11", "../data/worstcase-110-terms-4-m//inner-replacement54_11_21.dat", FULL);
relation* rel_reachable_5_1_2_3_4_5 = new relation(5, true, 5, 280, "rel_reachable_5_1_2_3_4_5", "../data/worstcase-110-terms-4-m//reachable_5_20.dat", FULL);
relation* rel_lam_2_0 = new relation(1, false, 2, 263, "rel_lam_2_0", "../data/worstcase-110-terms-4-m//lam_2_19.dat", FULL);
relation* rel_t_ret_to_var_10_10_9_8_7_6 = new relation(5, false, 10, 272, "rel_t_ret_to_var_10_10_9_8_7_6", "../data/worstcase-110-terms-4-m//t-ret-to-var_10_18.dat", FULL);
relation* rel_t_ret_to_ret_10_6_10_9_8_7 = new relation(5, false, 10, 285, "rel_t_ret_to_ret_10_6_10_9_8_7", "../data/worstcase-110-terms-4-m//t-ret-to-ret_10_17.dat", FULL);
relation* rel_free_2_1_2 = new relation(2, true, 2, 268, "rel_free_2_1_2", "../data/worstcase-110-terms-4-m//free_2_16.dat", FULL);
relation* rel_t_ret_to_ret_10_8_7_6 = new relation(3, false, 10, 285, "rel_t_ret_to_ret_10_8_7_6", "../data/worstcase-110-terms-4-m//t-ret-to-ret_10_15.dat", FULL);
relation* rel_free_2_2 = new relation(1, false, 2, 268, "rel_free_2_2", "../data/worstcase-110-terms-4-m//free_2_14.dat", FULL);
relation* rel_t_ret_to_var_10_1_2_3_4_5_6_7_8_9_10 = new relation(10, true, 10, 272, "rel_t_ret_to_var_10_1_2_3_4_5_6_7_8_9_10", "../data/worstcase-110-terms-4-m//t-ret-to-var_10_13.dat", FULL);
relation* rel_app_2_0 = new relation(1, false, 2, 259, "rel_app_2_0", "../data/worstcase-110-terms-4-m//app_2_12.dat", FULL);
relation* rel_inter_head1434_6_1_2_3_4_5_6 = new relation(6, true, 6, 271, "rel_inter_head1434_6_1_2_3_4_5_6", "../data/worstcase-110-terms-4-m//inter-head1434_6_11.dat", FULL);
relation* rel_inter_body1420_9_6 = new relation(1, false, 9, 275, "rel_inter_body1420_9_6", "../data/worstcase-110-terms-4-m//inter-body1420_9_10.dat", FULL);
relation* rel_inter_body1416_3_1_2 = new relation(2, false, 3, 284, "rel_inter_body1416_3_1_2", "../data/worstcase-110-terms-4-m//inter-body1416_3_9.dat", FULL);
relation* rel_t_ret_to_ret_10_1_2_3_4_5_6_7_8_9_10 = new relation(10, true, 10, 285, "rel_t_ret_to_ret_10_1_2_3_4_5_6_7_8_9_10", "../data/worstcase-110-terms-4-m//t-ret-to-ret_10_8.dat", FULL);
relation* rel_inter_body1416_3_1_2_3 = new relation(3, true, 3, 284, "rel_inter_body1416_3_1_2_3", "../data/worstcase-110-terms-4-m//inter-body1416_3_7.dat", FULL);
relation* rel_inter_head1428_7_1_2_3_4_5_6_7 = new relation(7, true, 7, 260, "rel_inter_head1428_7_1_2_3_4_5_6_7", "../data/worstcase-110-terms-4-m//inter-head1428_7_6.dat", FULL);
relation* rel_inter_body1420_9_1_2_3_4_5_6_7_8_9 = new relation(9, true, 9, 275, "rel_inter_body1420_9_1_2_3_4_5_6_7_8_9", "../data/worstcase-110-terms-4-m//inter-body1420_9_5.dat", FULL);
relation* rel_let_3_1_2_3 = new relation(3, true, 3, 282, "rel_let_3_1_2_3", "../data/worstcase-110-terms-4-m//let_3_4.dat", FULL);
relation* rel_let_3_0 = new relation(1, false, 3, 282, "rel_let_3_0", "../data/worstcase-110-terms-4-m//let_3_3.dat", FULL);
relation* rel_inter_head1437_13_1_2_3_4_5_6_7_8_9_10_11_12_13 = new relation(13, true, 13, 265, "rel_inter_head1437_13_1_2_3_4_5_6_7_8_9_10_11_12_13", "../data/worstcase-110-terms-4-m//inter-head1437_13_2.dat", FULL);
relation* rel_seq_2_2 = new relation(1, false, 2, 258, "rel_seq_2_2", "../data/worstcase-110-terms-4-m//seq_2_1.dat", FULL);
relation* rel_t_ret_to_ret_10_10_9_8_7_6 = new relation(5, false, 10, 285, "rel_t_ret_to_ret_10_10_9_8_7_6", "../data/worstcase-110-terms-4-m//t-ret-to-ret_10_0.dat", FULL);

RAM* scc4991 = new RAM(false, 1);
scc4991->add_relation(rel_ref_1_1, true);
scc4991->add_relation(rel_ref_1_0, true);
scc4991->add_rule(new parallel_acopy(rel_ref_1_0, rel_ref_1_1, DELTA, {1, 0}));

RAM* scc4992 = new RAM(false, 5);
scc4992->add_relation(rel_app_2_0, true);
scc4992->add_relation(rel_app_2_1_2, true);
scc4992->add_rule(new parallel_acopy(rel_app_2_0, rel_app_2_1_2, DELTA, {2, 0, 1}));

RAM* scc4993 = new RAM(true, 9);
scc4993->add_relation(rel_t_ret_to_ret_10_10_9_8_7_6, true);
scc4993->add_relation(rel_inter_head1437_13_1_2_3_4_5_6_7_8_9_10_11_12_13, true);
scc4993->add_relation(rel_let_3_0, false);
scc4993->add_relation(rel_inter_body1420_9_1_2_3_4_5_6_7_8_9, true);
scc4993->add_relation(rel_inter_head1428_7_1_2_3_4_5_6_7, true);
scc4993->add_relation(rel_t_ret_to_ret_10_1_2_3_4_5_6_7_8_9_10, true);
scc4993->add_relation(rel_inter_body1420_9_6, true);
scc4993->add_relation(rel_inter_head1434_6_1_2_3_4_5_6, true);
scc4993->add_relation(rel_app_2_0, false);
scc4993->add_relation(rel_t_ret_to_var_10_1_2_3_4_5_6_7_8_9_10, true);
scc4993->add_relation(rel_free_2_2, false);
scc4993->add_relation(rel_t_ret_to_ret_10_8_7_6, true);
scc4993->add_relation(rel_t_ret_to_ret_10_6_10_9_8_7, true);
scc4993->add_relation(rel_t_ret_to_var_10_10_9_8_7_6, true);
scc4993->add_relation(rel_lam_2_0, false);
scc4993->add_relation(rel_reachable_5_1_2_3_4_5, true);
scc4993->add_relation(rel_inner_replacement54_11_1_2_3_4_5_6_7_8_9_10_11, true);
scc4993->add_relation(rel_inner_replacement55_11_1_2_3_4_5_6_7_8_9_10_11, true);
scc4993->add_relation(rel_free_var_prop_9_1_2_3_4_5_6_7_8_9, true);
scc4993->add_relation(rel_ret_to_ret_10_1_2_3_4_5_6_7_8_9_10, true);
scc4993->add_relation(rel_var_to_ret_10_1_2_3_4_5_6_7_8_9_10, true);
scc4993->add_relation(rel_inner_replacement53_7_6_5_4_3_2, true);
scc4993->add_relation(rel_producer_5_5_4_3_2_1, true);
scc4993->add_relation(rel_ret_to_var_10_5_4_3_2_1, true);
scc4993->add_relation(rel_ret_to_ret_8_3_2_1, false);
scc4993->add_relation(rel_const_1_0, false);
scc4993->add_relation(rel_var_to_var_10_5_4_3_2_1, true);
scc4993->add_relation(rel_ret_to_ret_10_5_4_3_2_1, true);
scc4993->add_relation(rel_inner_replacement54_11_5_4_3_2_1, true);
scc4993->add_relation(rel_free_var_prop_9_1, true);
scc4993->add_relation(rel_var_to_var_10_1_2_3_4_5_6_7_8_9_10, true);
scc4993->add_relation(rel_inner_replacement55_11_1, true);
scc4993->add_relation(rel_inner_replacement53_7_1_2_3_4_5_6_7, true);
scc4993->add_relation(rel_var_to_ret_10_5_4_3_2_1, true);
scc4993->add_relation(rel_reachable_5_1, true);
scc4993->add_relation(rel_inter_head1431_6_1_2_3_4_5_6, true);
scc4993->add_relation(rel_ref_1_0, false);
scc4993->add_relation(rel_ret_to_var_10_1_2_3_4_5_6_7_8_9_10, true);
scc4993->add_relation(rel_seq_2_0, false);
scc4993->add_rule(new parallel_acopy(rel_reachable_5_1, rel_reachable_5_1_2_3_4_5, DELTA, {0, 5, 1, 2, 3, 4}));
scc4993->add_rule(new parallel_acopy(rel_inner_replacement53_7_6_5_4_3_2, rel_inner_replacement53_7_1_2_3_4_5_6_7, DELTA, {5, 4, 3, 2, 1, 7, 0, 6}));
scc4993->add_rule(new parallel_join(rel_inner_replacement54_11_1_2_3_4_5_6_7_8_9_10_11, rel_inner_replacement53_7_6_5_4_3_2, DELTA, rel_t_ret_to_ret_10_6_10_9_8_7, DELTA, {9, 10, 11, 12, 13, 4, 3, 2, 1, 6, 7}));
scc4993->add_rule(new parallel_join(rel_inner_replacement53_7_1_2_3_4_5_6_7, rel_app_2_0, FULL, rel_reachable_5_1, DELTA, {0, 4, 5, 6, 7, 1, 2}));
scc4993->add_rule(new parallel_join(rel_t_ret_to_var_10_1_2_3_4_5_6_7_8_9_10, rel_ret_to_var_10_5_4_3_2_1, DELTA, rel_t_ret_to_ret_10_10_9_8_7_6, FULL, {12, 13, 14, 15, 16, 6, 7, 8, 9, 10}));
scc4993->add_rule(new parallel_acopy(rel_ret_to_ret_10_5_4_3_2_1, rel_ret_to_ret_10_1_2_3_4_5_6_7_8_9_10, DELTA, {4, 3, 2, 1, 0, 10, 5, 6, 7, 8, 9}));
scc4993->add_rule(new parallel_join(rel_t_ret_to_var_10_1_2_3_4_5_6_7_8_9_10, rel_ret_to_var_10_5_4_3_2_1, DELTA, rel_t_ret_to_ret_10_10_9_8_7_6, DELTA, {12, 13, 14, 15, 16, 6, 7, 8, 9, 10}));
scc4993->add_rule(new parallel_join(rel_inter_head1431_6_1_2_3_4_5_6, rel_seq_2_0, FULL, rel_reachable_5_1, DELTA, {1, 2, 6, 5, 7, 4}));
scc4993->add_rule(new parallel_acopy(rel_inner_replacement54_11_5_4_3_2_1, rel_inner_replacement54_11_1_2_3_4_5_6_7_8_9_10_11, DELTA, {4, 3, 2, 1, 0, 11, 5, 6, 7, 8, 9, 10}));
scc4993->add_rule(new parallel_copy(rel_reachable_5_1_2_3_4_5, rel_inter_head1434_6_1_2_3_4_5_6, DELTA, {0, 5, 3, 1, 4}));
scc4993->add_rule(new parallel_copy(rel_reachable_5_1_2_3_4_5, rel_inter_head1431_6_1_2_3_4_5_6, DELTA, {0, 5, 3, 2, 4}));
scc4993->add_rule(new parallel_join(rel_t_ret_to_var_10_1_2_3_4_5_6_7_8_9_10, rel_producer_5_5_4_3_2_1, DELTA, rel_ret_to_var_10_5_4_3_2_1, DELTA, {4, 3, 2, 1, 0, 7, 8, 9, 10, 11}));
scc4993->add_rule(new parallel_join(rel_t_ret_to_var_10_1_2_3_4_5_6_7_8_9_10, rel_var_to_var_10_5_4_3_2_1, FULL, rel_t_ret_to_var_10_10_9_8_7_6, DELTA, {12, 13, 14, 15, 16, 6, 7, 8, 9, 10}));
scc4993->add_rule(new parallel_join(rel_t_ret_to_var_10_1_2_3_4_5_6_7_8_9_10, rel_var_to_var_10_5_4_3_2_1, DELTA, rel_t_ret_to_var_10_10_9_8_7_6, DELTA, {12, 13, 14, 15, 16, 6, 7, 8, 9, 10}));
scc4993->add_rule(new parallel_join(rel_producer_5_5_4_3_2_1, rel_lam_2_0, FULL, rel_reachable_5_1, DELTA, {7, 6, 5, 4, 0}));
scc4993->add_rule(new parallel_join(rel_t_ret_to_ret_10_1_2_3_4_5_6_7_8_9_10, rel_producer_5_5_4_3_2_1, DELTA, rel_ret_to_ret_10_5_4_3_2_1, DELTA, {4, 3, 2, 1, 0, 7, 8, 9, 10, 11}));
scc4993->add_rule(new parallel_acopy(rel_var_to_ret_10_5_4_3_2_1, rel_var_to_ret_10_1_2_3_4_5_6_7_8_9_10, DELTA, {4, 3, 2, 1, 0, 10, 5, 6, 7, 8, 9}));
scc4993->add_rule(new parallel_join(rel_inter_head1434_6_1_2_3_4_5_6, rel_app_2_0, FULL, rel_reachable_5_1, DELTA, {2, 6, 1, 5, 7, 4}));
scc4993->add_rule(new parallel_join(rel_t_ret_to_var_10_1_2_3_4_5_6_7_8_9_10, rel_var_to_var_10_5_4_3_2_1, DELTA, rel_t_ret_to_var_10_10_9_8_7_6, FULL, {12, 13, 14, 15, 16, 6, 7, 8, 9, 10}));
scc4993->add_rule(new parallel_join(rel_var_to_var_10_1_2_3_4_5_6_7_8_9_10, rel_free_2_2, FULL, rel_inter_body1420_9_6, DELTA, {2, 6, 5, 8, 7, 2, 11, 9, 4, 10}));
scc4993->add_rule(new parallel_acopy(rel_t_ret_to_ret_10_8_7_6, rel_t_ret_to_ret_10_1_2_3_4_5_6_7_8_9_10, DELTA, {7, 6, 5, 10, 0, 1, 2, 3, 4, 8, 9}));
scc4993->add_rule(new parallel_copy(rel_reachable_5_1_2_3_4_5, rel_inter_head1437_13_1_2_3_4_5_6_7_8_9_10_11_12_13, DELTA, {12, 10, 11, 8, 2}));
scc4993->add_rule(new parallel_acopy(rel_t_ret_to_ret_10_6_10_9_8_7, rel_t_ret_to_ret_10_1_2_3_4_5_6_7_8_9_10, DELTA, {5, 9, 8, 7, 6, 10, 0, 1, 2, 3, 4}));
scc4993->add_rule(new parallel_copy(rel_ret_to_ret_10_1_2_3_4_5_6_7_8_9_10, rel_producer_5_5_4_3_2_1, DELTA, {4, 3, 2, 1, 0, 4, 3, 2, 1, 0}));
scc4993->add_rule(new parallel_acopy(rel_t_ret_to_var_10_10_9_8_7_6, rel_t_ret_to_var_10_1_2_3_4_5_6_7_8_9_10, DELTA, {9, 8, 7, 6, 5, 10, 0, 1, 2, 3, 4}));
scc4993->add_rule(new parallel_join(rel_inner_replacement55_11_1_2_3_4_5_6_7_8_9_10_11, rel_producer_5_5_4_3_2_1, FULL, rel_inner_replacement54_11_5_4_3_2_1, DELTA, {4, 3, 2, 1, 0, 7, 8, 9, 10, 11, 12}));
scc4993->add_rule(new parallel_join(rel_inner_replacement55_11_1_2_3_4_5_6_7_8_9_10_11, rel_producer_5_5_4_3_2_1, DELTA, rel_inner_replacement54_11_5_4_3_2_1, DELTA, {4, 3, 2, 1, 0, 7, 8, 9, 10, 11, 12}));
scc4993->add_rule(new parallel_join(rel_inner_replacement54_11_1_2_3_4_5_6_7_8_9_10_11, rel_inner_replacement53_7_6_5_4_3_2, FULL, rel_t_ret_to_ret_10_6_10_9_8_7, DELTA, {9, 10, 11, 12, 13, 4, 3, 2, 1, 6, 7}));
scc4993->add_rule(new parallel_copy(rel_ret_to_ret_10_1_2_3_4_5_6_7_8_9_10, rel_inter_head1437_13_1_2_3_4_5_6_7_8_9_10_11_12_13, DELTA, {12, 10, 11, 8, 2, 10, 11, 8, 2, 9}));
scc4993->add_rule(new parallel_join(rel_t_ret_to_ret_10_1_2_3_4_5_6_7_8_9_10, rel_t_ret_to_var_10_10_9_8_7_6, DELTA, rel_var_to_ret_10_5_4_3_2_1, DELTA, {6, 7, 8, 9, 10, 12, 13, 14, 15, 16}));
scc4993->add_rule(new parallel_join(rel_t_ret_to_ret_10_1_2_3_4_5_6_7_8_9_10, rel_t_ret_to_var_10_10_9_8_7_6, FULL, rel_var_to_ret_10_5_4_3_2_1, DELTA, {6, 7, 8, 9, 10, 12, 13, 14, 15, 16}));
scc4993->add_rule(new parallel_join(rel_inter_head1428_7_1_2_3_4_5_6_7, rel_let_3_0, FULL, rel_reachable_5_1, DELTA, {7, 2, 1, 6, 8, 5, 3}));
scc4993->add_rule(new parallel_copy(rel_ret_to_var_10_1_2_3_4_5_6_7_8_9_10, rel_inter_head1428_7_1_2_3_4_5_6_7, DELTA, {1, 5, 3, 0, 4, 2, 5, 3, 0, 4}));
scc4993->add_rule(new parallel_copy(rel_reachable_5_1_2_3_4_5, rel_inter_head1428_7_1_2_3_4_5_6_7, DELTA, {1, 5, 3, 0, 4}));
scc4993->add_rule(new parallel_join(rel_t_ret_to_var_10_1_2_3_4_5_6_7_8_9_10, rel_ret_to_var_10_5_4_3_2_1, FULL, rel_t_ret_to_ret_10_10_9_8_7_6, DELTA, {12, 13, 14, 15, 16, 6, 7, 8, 9, 10}));
scc4993->add_rule(new parallel_copy(rel_free_var_prop_9_1_2_3_4_5_6_7_8_9, rel_inter_head1437_13_1_2_3_4_5_6_7_8_9_10_11_12_13, DELTA, {1, 4, 3, 6, 5, 10, 11, 8, 2}));
scc4993->add_rule(new parallel_join(rel_inner_replacement55_11_1_2_3_4_5_6_7_8_9_10_11, rel_producer_5_5_4_3_2_1, DELTA, rel_inner_replacement54_11_5_4_3_2_1, FULL, {4, 3, 2, 1, 0, 7, 8, 9, 10, 11, 12}));
scc4993->add_rule(new parallel_acopy(rel_t_ret_to_ret_10_10_9_8_7_6, rel_t_ret_to_ret_10_1_2_3_4_5_6_7_8_9_10, DELTA, {9, 8, 7, 6, 5, 10, 0, 1, 2, 3, 4}));
scc4993->add_rule(new parallel_join(rel_t_ret_to_ret_10_1_2_3_4_5_6_7_8_9_10, rel_ret_to_ret_8_3_2_1, FULL, rel_t_ret_to_ret_10_8_7_6, DELTA, {10, 11, 12, 13, 14, 4, 5, 6, 7, 8}));
scc4993->add_rule(new parallel_acopy(rel_inter_body1420_9_6, rel_inter_body1420_9_1_2_3_4_5_6_7_8_9, DELTA, {5, 9, 0, 1, 2, 3, 4, 6, 7, 8}));
scc4993->add_rule(new parallel_join(rel_t_ret_to_var_10_1_2_3_4_5_6_7_8_9_10, rel_producer_5_5_4_3_2_1, FULL, rel_ret_to_var_10_5_4_3_2_1, DELTA, {4, 3, 2, 1, 0, 7, 8, 9, 10, 11}));
scc4993->add_rule(new parallel_copy(rel_reachable_5_1_2_3_4_5, rel_inter_head1434_6_1_2_3_4_5_6, DELTA, {2, 5, 3, 1, 4}));
scc4993->add_rule(new parallel_join(rel_inter_head1437_13_1_2_3_4_5_6_7_8_9_10_11_12_13, rel_lam_2_0, FULL, rel_inner_replacement55_11_1, DELTA, {13, 0, 10, 5, 4, 7, 6, 1, 9, 11, 12, 8, 2}));
scc4993->add_rule(new parallel_acopy(rel_ret_to_var_10_5_4_3_2_1, rel_ret_to_var_10_1_2_3_4_5_6_7_8_9_10, DELTA, {4, 3, 2, 1, 0, 10, 5, 6, 7, 8, 9}));
scc4993->add_rule(new parallel_copy(rel_ret_to_var_10_1_2_3_4_5_6_7_8_9_10, rel_inter_head1437_13_1_2_3_4_5_6_7_8_9_10_11_12_13, DELTA, {0, 11, 8, 2, 9, 7, 10, 11, 8, 2}));
scc4993->add_rule(new parallel_join(rel_t_ret_to_ret_10_1_2_3_4_5_6_7_8_9_10, rel_producer_5_5_4_3_2_1, DELTA, rel_ret_to_ret_10_5_4_3_2_1, FULL, {4, 3, 2, 1, 0, 7, 8, 9, 10, 11}));
scc4993->add_rule(new parallel_join(rel_producer_5_5_4_3_2_1, rel_const_1_0, FULL, rel_reachable_5_1, DELTA, {6, 5, 4, 3, 0}));
scc4993->add_rule(new parallel_join(rel_inner_replacement54_11_1_2_3_4_5_6_7_8_9_10_11, rel_inner_replacement53_7_6_5_4_3_2, DELTA, rel_t_ret_to_ret_10_6_10_9_8_7, FULL, {9, 10, 11, 12, 13, 4, 3, 2, 1, 6, 7}));
scc4993->add_rule(new parallel_acopy(rel_inner_replacement55_11_1, rel_inner_replacement55_11_1_2_3_4_5_6_7_8_9_10_11, DELTA, {0, 11, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}));
scc4993->add_rule(new parallel_join(rel_inter_body1420_9_1_2_3_4_5_6_7_8_9, rel_lam_2_0, FULL, rel_free_var_prop_9_1, DELTA, {10, 5, 4, 7, 6, 0, 9, 11, 8}));
scc4993->add_rule(new parallel_acopy(rel_var_to_var_10_5_4_3_2_1, rel_var_to_var_10_1_2_3_4_5_6_7_8_9_10, DELTA, {4, 3, 2, 1, 0, 10, 5, 6, 7, 8, 9}));
scc4993->add_rule(new parallel_join(rel_var_to_ret_10_1_2_3_4_5_6_7_8_9_10, rel_ref_1_0, FULL, rel_reachable_5_1, DELTA, {1, 3, 4, 5, 6, 0, 3, 4, 5, 6}));
scc4993->add_rule(new parallel_copy(rel_reachable_5_1_2_3_4_5, rel_inter_head1431_6_1_2_3_4_5_6, DELTA, {1, 5, 3, 2, 4}));
scc4993->add_rule(new parallel_join(rel_t_ret_to_ret_10_1_2_3_4_5_6_7_8_9_10, rel_t_ret_to_var_10_10_9_8_7_6, DELTA, rel_var_to_ret_10_5_4_3_2_1, FULL, {6, 7, 8, 9, 10, 12, 13, 14, 15, 16}));
scc4993->add_rule(new parallel_copy(rel_reachable_5_1_2_3_4_5, rel_inter_head1428_7_1_2_3_4_5_6_7, DELTA, {6, 5, 3, 0, 4}));
scc4993->add_rule(new parallel_join(rel_t_ret_to_var_10_1_2_3_4_5_6_7_8_9_10, rel_producer_5_5_4_3_2_1, DELTA, rel_ret_to_var_10_5_4_3_2_1, FULL, {4, 3, 2, 1, 0, 7, 8, 9, 10, 11}));
scc4993->add_rule(new parallel_acopy(rel_free_var_prop_9_1, rel_free_var_prop_9_1_2_3_4_5_6_7_8_9, DELTA, {0, 9, 1, 2, 3, 4, 5, 6, 7, 8}));
scc4993->add_rule(new parallel_join(rel_t_ret_to_ret_10_1_2_3_4_5_6_7_8_9_10, rel_producer_5_5_4_3_2_1, FULL, rel_ret_to_ret_10_5_4_3_2_1, DELTA, {4, 3, 2, 1, 0, 7, 8, 9, 10, 11}));

RAM* scc4994 = new RAM(false, 13);
scc4994->add_relation(rel_seq_2_2, true);
scc4994->add_relation(rel_seq_2_1_2, true);
scc4994->add_rule(new parallel_acopy(rel_seq_2_2, rel_seq_2_1_2, DELTA, {1, 2, 0}));

RAM* scc4995 = new RAM(false, 18);
scc4995->add_relation(rel_let_3_0, true);
scc4995->add_relation(rel_let_3_1_2_3, true);
scc4995->add_rule(new parallel_acopy(rel_let_3_0, rel_let_3_1_2_3, DELTA, {3, 0, 1, 2}));

RAM* scc4996 = new RAM(false, 3);
scc4996->add_relation(rel_const_1_0, true);
scc4996->add_relation(rel_const_1_1, true);
scc4996->add_rule(new parallel_acopy(rel_const_1_0, rel_const_1_1, DELTA, {1, 0}));

RAM* scc4997 = new RAM(false, 7);
scc4997->add_relation(rel_reachable_5_1_2_3_4_5, true);
scc4997->add_relation(rel_program_1_1, false);
scc4997->add_rule(new parallel_copy(rel_reachable_5_1_2_3_4_5, rel_program_1_1, FULL, {0, 0, 0, 0, 0}));

RAM* scc4998 = new RAM(false, 11);
scc4998->add_relation(rel_ret_to_ret_8_1_2_3_4_5_6_7_8, true);
scc4998->add_relation(rel_ret_to_ret_8_3_2_1, true);
scc4998->add_rule(new parallel_acopy(rel_ret_to_ret_8_3_2_1, rel_ret_to_ret_8_1_2_3_4_5_6_7_8, DELTA, {2, 1, 0, 8, 3, 4, 5, 6, 7}));

RAM* scc4999 = new RAM(false, 15);
scc4999->add_relation(rel_lam_2_1_2, true);
scc4999->add_relation(rel_lam_2_2, true);
scc4999->add_rule(new parallel_acopy(rel_lam_2_2, rel_lam_2_1_2, DELTA, {1, 2, 0}));

RAM* scc5000 = new RAM(false, 17);
scc5000->add_relation(rel_seq_2_1_2, true);
scc5000->add_relation(rel_seq_2_1, true);
scc5000->add_rule(new parallel_acopy(rel_seq_2_1, rel_seq_2_1_2, DELTA, {0, 2, 1}));

RAM* scc5001 = new RAM(false, 2);
scc5001->add_relation(rel_lam_2_0, true);
scc5001->add_relation(rel_lam_2_1_2, true);
scc5001->add_rule(new parallel_acopy(rel_lam_2_0, rel_lam_2_1_2, DELTA, {2, 0, 1}));

RAM* scc5002 = new RAM(false, 6);
scc5002->add_relation(rel_let_3_1_2_3, true);
scc5002->add_relation(rel_let_3_3, true);
scc5002->add_rule(new parallel_acopy(rel_let_3_3, rel_let_3_1_2_3, DELTA, {2, 3, 0, 1}));

RAM* scc5003 = new RAM(true, 10);
scc5003->add_relation(rel_seq_2_2, false);
scc5003->add_relation(rel_inter_body1416_3_1_2_3, true);
scc5003->add_relation(rel_inter_body1416_3_1_2, true);
scc5003->add_relation(rel_free_2_2, true);
scc5003->add_relation(rel_free_2_1_2, true);
scc5003->add_relation(rel_let_3_2, false);
scc5003->add_relation(rel_app_2_2, false);
scc5003->add_relation(rel_lam_2_2, false);
scc5003->add_relation(rel_seq_2_1, false);
scc5003->add_relation(rel_inter_body1424_3_1_2_3, true);
scc5003->add_relation(rel_app_2_1, false);
scc5003->add_relation(rel_let_3_3, false);
scc5003->add_relation(rel_inter_body1424_3_2_3, true);
scc5003->add_rule(new parallel_join(rel_free_2_1_2, rel_app_2_2, FULL, rel_free_2_2, DELTA, {4, 1}));
scc5003->add_rule(new parallel_acopy(rel_free_2_2, rel_free_2_1_2, DELTA, {1, 2, 0}));
scc5003->add_rule(new parallel_copy_filter(rel_free_2_1_2, rel_inter_body1416_3_1_2, DELTA, {1, 3}, [](const u64* const data){ return !(data[0] == data[1]); }));
scc5003->add_rule(new parallel_join(rel_free_2_1_2, rel_free_2_2, DELTA, rel_let_3_2, FULL, {2, 3}));
scc5003->add_rule(new parallel_acopy(rel_inter_body1424_3_2_3, rel_inter_body1424_3_1_2_3, DELTA, {1, 2, 3, 0}));
scc5003->add_rule(new parallel_join(rel_inter_body1424_3_1_2_3, rel_free_2_2, DELTA, rel_let_3_3, FULL, {3, 4, 2}));
scc5003->add_rule(new parallel_copy_filter(rel_free_2_1_2, rel_inter_body1424_3_2_3, DELTA, {1, 3}, [](const u64* const data){ return !(data[0] == data[1]); }));
scc5003->add_rule(new parallel_join(rel_free_2_1_2, rel_app_2_1, FULL, rel_free_2_2, DELTA, {4, 1}));
scc5003->add_rule(new parallel_acopy(rel_inter_body1416_3_1_2, rel_inter_body1416_3_1_2_3, DELTA, {0, 1, 3, 2}));
scc5003->add_rule(new parallel_join(rel_free_2_1_2, rel_free_2_2, DELTA, rel_seq_2_2, FULL, {2, 3}));
scc5003->add_rule(new parallel_join(rel_inter_body1416_3_1_2_3, rel_lam_2_2, FULL, rel_free_2_2, DELTA, {2, 4, 1}));
scc5003->add_rule(new parallel_join(rel_free_2_1_2, rel_free_2_2, DELTA, rel_seq_2_1, FULL, {2, 3}));

RAM* scc5004 = new RAM(false, 14);
scc5004->add_relation(rel_free_2_1_2, true);
scc5004->add_relation(rel_ref_1_1, false);
scc5004->add_rule(new parallel_copy(rel_free_2_1_2, rel_ref_1_1, FULL, {0, 1}));

RAM* scc5005 = new RAM(false, 19);
scc5005->add_relation(rel_app_2_1, true);
scc5005->add_relation(rel_app_2_1_2, true);
scc5005->add_rule(new parallel_acopy(rel_app_2_1, rel_app_2_1_2, DELTA, {0, 2, 1}));

RAM* scc5006 = new RAM(false, 4);
scc5006->add_relation(rel_app_2_2, true);
scc5006->add_relation(rel_app_2_1_2, true);
scc5006->add_rule(new parallel_acopy(rel_app_2_2, rel_app_2_1_2, DELTA, {1, 2, 0}));

RAM* scc5007 = new RAM(false, 8);
scc5007->add_relation(rel_seq_2_1_2, true);
scc5007->add_relation(rel_seq_2_0, true);
scc5007->add_rule(new parallel_acopy(rel_seq_2_0, rel_seq_2_1_2, DELTA, {2, 0, 1}));

RAM* scc5008 = new RAM(false, 12);
scc5008->add_relation(rel_let_3_1_2_3, true);
scc5008->add_relation(rel_let_3_2, true);
scc5008->add_rule(new parallel_acopy(rel_let_3_2, rel_let_3_1_2_3, DELTA, {1, 3, 0, 2}));

RAM* scc5009 = new RAM(false, 16);
scc5009->add_relation(rel_inter_head1437_13_1_2_3_4_5_6_7_8_9_10_11_12_13, false);
scc5009->add_relation(rel_app_step_10_1_2_3_4_5_6_7_8_9_10, true);
scc5009->add_rule(new parallel_copy(rel_app_step_10_1_2_3_4_5_6_7_8_9_10, rel_inter_head1437_13_1_2_3_4_5_6_7_8_9_10_11_12_13, FULL, {10, 11, 8, 2, 9, 12, 10, 11, 8, 2}));

LIE* lie = new LIE();
lie->add_relation(rel_seq_2_0);
lie->add_relation(rel_ret_to_var_10_1_2_3_4_5_6_7_8_9_10);
lie->add_relation(rel_app_2_1_2);
lie->add_relation(rel_inter_body1424_3_2_3);
lie->add_relation(rel_let_3_3);
lie->add_relation(rel_ref_1_0);
lie->add_relation(rel_inter_head1431_6_1_2_3_4_5_6);
lie->add_relation(rel_reachable_5_1);
lie->add_relation(rel_var_to_ret_10_5_4_3_2_1);
lie->add_relation(rel_const_1_1);
lie->add_relation(rel_inner_replacement53_7_1_2_3_4_5_6_7);
lie->add_relation(rel_inner_replacement55_11_1);
lie->add_relation(rel_var_to_var_10_1_2_3_4_5_6_7_8_9_10);
lie->add_relation(rel_app_2_1);
lie->add_relation(rel_free_var_prop_9_1);
lie->add_relation(rel_inter_body1424_3_1_2_3);
lie->add_relation(rel_seq_2_1);
lie->add_relation(rel_lam_2_2);
lie->add_relation(rel_seq_2_1_2);
lie->add_relation(rel_program_1_1);
lie->add_relation(rel_inner_replacement54_11_5_4_3_2_1);
lie->add_relation(rel_ret_to_ret_10_5_4_3_2_1);
lie->add_relation(rel_var_to_var_10_5_4_3_2_1);
lie->add_relation(rel_const_1_0);
lie->add_relation(rel_ret_to_ret_8_3_2_1);
lie->add_relation(rel_ret_to_var_10_5_4_3_2_1);
lie->add_relation(rel_producer_5_5_4_3_2_1);
lie->add_relation(rel_app_2_2);
lie->add_relation(rel_inner_replacement53_7_6_5_4_3_2);
lie->add_relation(rel_let_3_2);
lie->add_relation(rel_var_to_ret_10_1_2_3_4_5_6_7_8_9_10);
lie->add_relation(rel_ret_to_ret_10_1_2_3_4_5_6_7_8_9_10);
lie->add_relation(rel_ret_to_ret_8_1_2_3_4_5_6_7_8);
lie->add_relation(rel_free_var_prop_9_1_2_3_4_5_6_7_8_9);
lie->add_relation(rel_lam_2_1_2);
lie->add_relation(rel_inner_replacement55_11_1_2_3_4_5_6_7_8_9_10_11);
lie->add_relation(rel_app_step_10_1_2_3_4_5_6_7_8_9_10);
lie->add_relation(rel_ref_1_1);
lie->add_relation(rel_inner_replacement54_11_1_2_3_4_5_6_7_8_9_10_11);
lie->add_relation(rel_reachable_5_1_2_3_4_5);
lie->add_relation(rel_lam_2_0);
lie->add_relation(rel_t_ret_to_var_10_10_9_8_7_6);
lie->add_relation(rel_t_ret_to_ret_10_6_10_9_8_7);
lie->add_relation(rel_free_2_1_2);
lie->add_relation(rel_t_ret_to_ret_10_8_7_6);
lie->add_relation(rel_free_2_2);
lie->add_relation(rel_t_ret_to_var_10_1_2_3_4_5_6_7_8_9_10);
lie->add_relation(rel_app_2_0);
lie->add_relation(rel_inter_head1434_6_1_2_3_4_5_6);
lie->add_relation(rel_inter_body1420_9_6);
lie->add_relation(rel_inter_body1416_3_1_2);
lie->add_relation(rel_t_ret_to_ret_10_1_2_3_4_5_6_7_8_9_10);
lie->add_relation(rel_inter_body1416_3_1_2_3);
lie->add_relation(rel_inter_head1428_7_1_2_3_4_5_6_7);
lie->add_relation(rel_inter_body1420_9_1_2_3_4_5_6_7_8_9);
lie->add_relation(rel_let_3_1_2_3);
lie->add_relation(rel_let_3_0);
lie->add_relation(rel_inter_head1437_13_1_2_3_4_5_6_7_8_9_10_11_12_13);
lie->add_relation(rel_seq_2_2);
lie->add_relation(rel_t_ret_to_ret_10_10_9_8_7_6);
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
lie->add_scc_dependance(scc4991, scc4993);
lie->add_scc_dependance(scc4992, scc4993);
lie->add_scc_dependance(scc4993, scc5009);
lie->add_scc_dependance(scc4994, scc5003);
lie->add_scc_dependance(scc4995, scc4993);
lie->add_scc_dependance(scc4996, scc4993);
lie->add_scc_dependance(scc4997, scc4993);
lie->add_scc_dependance(scc4998, scc4993);
lie->add_scc_dependance(scc4999, scc5003);
lie->add_scc_dependance(scc5000, scc5003);
lie->add_scc_dependance(scc5001, scc4993);
lie->add_scc_dependance(scc5002, scc5003);
lie->add_scc_dependance(scc5003, scc4993);
lie->add_scc_dependance(scc5004, scc5003);
lie->add_scc_dependance(scc5005, scc5003);
lie->add_scc_dependance(scc5006, scc5003);
lie->add_scc_dependance(scc5007, scc4993);
lie->add_scc_dependance(scc5008, scc5003);



  // Enable IO
 //lie->enable_share_io();
 //lie->set_output_dir(../data/worstcase-110-terms-4-m-output); // Write to this directory
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
