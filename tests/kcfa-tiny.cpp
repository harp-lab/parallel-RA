#include "../src/parallel_RA_inc.h"


int main(int argc, char **argv)
{
    mpi_comm mcomm;
    mcomm.create(argc, argv);


    relation* rel_Lam_4_1_2_3_4 = new relation(4, true, 4, 261, "rel_Lam_4_1_2_3_4", "../data/g4470/Lam_4_38", FULL);
    relation* rel_AEval_14_7_6_5_4_3_2_1 = new relation(7, false, 14, 274, "rel_AEval_14_7_6_5_4_3_2_1", "../data/g4470/AEval_14_37", FULL);
    relation* rel_App_4_1_2_3_4 = new relation(4, true, 4, 271, "rel_App_4_1_2_3_4", "../data/g4470/App_4_36", FULL);
    relation* rel_INT1_20_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20 = new relation(20, true, 20, 260, "rel_INT1_20_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20", "../data/g4470/INT1_20_35", FULL);
    relation* rel_INT00_10_1_2_3_4_5_6_7_8_9_10 = new relation(10, true, 10, 265, "rel_INT00_10_1_2_3_4_5_6_7_8_9_10", "../data/g4470/INT00_10_34", FULL);
    relation* rel_INT_A_15_15_14_13_12_11_10 = new relation(6, false, 15, 258, "rel_INT_A_15_15_14_13_12_11_10", "../data/g4470/INT_A_15_33", FULL);
    relation* rel_INT2_27_10_9_8_7_6_5_3 = new relation(7, false, 27, 268, "rel_INT2_27_10_9_8_7_6_5_3", "../data/g4470/INT2_27_32", FULL);
    relation* rel_Lam_4_ = new relation(0, false, 4, 261, "rel_Lam_4_", "../data/g4470/Lam_4_31", FULL);
    relation* rel_INT_D_16_10_2 = new relation(2, false, 16, 264, "rel_INT_D_16_10_2", "../data/g4470/INT_D_16_30", FULL);
    relation* rel_Prog_1_ = new relation(0, false, 1, 259, "rel_Prog_1_", "../data/g4470/Prog_1_29", FULL);
    relation* rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31 = new relation(31, true, 31, 267, "rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31", "../data/g4470/APP_31_28", FULL);
    relation* rel_INT0_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17 = new relation(17, true, 17, 256, "rel_INT0_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17", "../data/g4470/INT0_17_27", FULL);
    relation* rel_INT_C_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17 = new relation(17, true, 17, 272, "rel_INT_C_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17", "../data/g4470/INT_C_17_26", FULL);
    relation* rel_INT_C_17_2_10 = new relation(2, false, 17, 272, "rel_INT_C_17_2_10", "../data/g4470/INT_C_17_25", FULL);
    relation* rel_Time_6_ = new relation(0, false, 6, 273, "rel_Time_6_", "../data/g4470/Time_6_24", FULL);
    relation* rel_INT_D_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16 = new relation(16, true, 16, 264, "rel_INT_D_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16", "../data/g4470/INT_D_16_23", FULL);
    relation* rel_INT1_20_10_9_8_7_6_5_4 = new relation(7, false, 20, 260, "rel_INT1_20_10_9_8_7_6_5_4", "../data/g4470/INT1_20_22", FULL);
    relation* rel_INT00_10_10_9_8_7_6_5_2 = new relation(7, false, 10, 265, "rel_INT00_10_10_9_8_7_6_5_2", "../data/g4470/INT00_10_21", FULL);
    relation* rel_AEval_14_7_6_5_4_3_2 = new relation(6, false, 14, 274, "rel_AEval_14_7_6_5_4_3_2", "../data/g4470/AEval_14_20", FULL);
    relation* rel_ReachesClo_7_1_2_3_4_5_6_7 = new relation(7, true, 7, 263, "rel_ReachesClo_7_1_2_3_4_5_6_7", "../data/g4470/ReachesClo_7_19", FULL);
    relation* rel_Time_6_1_2_3_4_5_6 = new relation(6, true, 6, 273, "rel_Time_6_1_2_3_4_5_6", "../data/g4470/Time_6_18", FULL);
    relation* rel_ReachesCfg_7_1 = new relation(1, false, 7, 257, "rel_ReachesCfg_7_1", "../data/g4470/ReachesCfg_7_17", FULL);
    relation* rel_Store_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14 = new relation(14, true, 14, 270, "rel_Store_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14", "../data/g4470/Store_14_16", FULL);
    relation* rel_INT2_27_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27 = new relation(27, true, 27, 268, "rel_INT2_27_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27", "../data/g4470/INT2_27_15", FULL);
    relation* rel_INT_A_15_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15 = new relation(15, true, 15, 258, "rel_INT_A_15_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15", "../data/g4470/INT_A_15_14", FULL);
    relation* rel_Lam_4_1 = new relation(1, false, 4, 261, "rel_Lam_4_1", "../data/g4470/Lam_4_13", FULL);
    relation* rel_Step_21_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21 = new relation(21, true, 21, 275, "rel_Step_21_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21", "../data/g4470/Step_21_12", FULL);
    relation* rel_Step_21_15 = new relation(1, false, 21, 275, "rel_Step_21_15", "../data/g4470/Step_21_11", FULL);
    relation* rel_Prog_1_1 = new relation(1, true, 1, 259, "rel_Prog_1_1", "../data/g4470/Prog_1_10", FULL);
    relation* rel_ReachesCfg_7_1_2_3_4_5_6_7 = new relation(7, true, 7, 257, "rel_ReachesCfg_7_1_2_3_4_5_6_7", "../data/g4470/ReachesCfg_7_9", FULL);
    relation* rel_Var_2_2 = new relation(1, false, 2, 266, "rel_Var_2_2", "../data/g4470/Var_2_8", FULL);
    relation* rel_INT0_17_11 = new relation(1, false, 17, 256, "rel_INT0_17_11", "../data/g4470/INT0_17_7", FULL);
    relation* rel_Var_2_1 = new relation(1, false, 2, 266, "rel_Var_2_1", "../data/g4470/Var_2_6", FULL);
    relation* rel_Var_2_1_2 = new relation(2, true, 2, 266, "rel_Var_2_1_2", "../data/g4470/Var_2_5", FULL);
    relation* rel_INT_B_16_1 = new relation(1, false, 16, 262, "rel_INT_B_16_1", "../data/g4470/INT_B_16_4", FULL);
    relation* rel_INT_B_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16 = new relation(16, true, 16, 262, "rel_INT_B_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16", "../data/g4470/INT_B_16_3", FULL);
    relation* rel_App_4_1 = new relation(1, false, 4, 271, "rel_App_4_1", "../data/g4470/App_4_2", FULL);
    relation* rel_Store_14_1 = new relation(1, false, 14, 270, "rel_Store_14_1", "../data/g4470/Store_14_1", FULL);
    relation* rel_AEval_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14 = new relation(14, true, 14, 274, "rel_AEval_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14", "../data/g4470/AEval_14_0", FULL);

    RAM* scc4471 = new RAM(false, 1);
    scc4471->add_relation(rel_Var_2_1_2, true);
    scc4471->add_relation(rel_Var_2_1, true);
    scc4471->add_rule(new parallel_acopy(rel_Var_2_1, rel_Var_2_1_2, DELTA, {0, 2, 1}));

    RAM* scc4472 = new RAM(false, 5);
    scc4472->add_relation(rel_ReachesClo_7_1_2_3_4_5_6_7, true);
    scc4472->add_relation(rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31, false);
    scc4472->add_rule(new parallel_copy(rel_ReachesClo_7_1_2_3_4_5_6_7, rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31, FULL, {7, 8, 9, 10, 11, 12, 13}));

    RAM* scc4473 = new RAM(false, 9);
    scc4473->add_relation(rel_Var_2_1_2, true);
    scc4473->add_relation(rel_Var_2_2, true);
    scc4473->add_rule(new parallel_acopy(rel_Var_2_2, rel_Var_2_1_2, DELTA, {1, 2, 0}));

    RAM* scc4474 = new RAM(false, 3);
    scc4474->add_relation(rel_Prog_1_1, true);
    scc4474->add_relation(rel_Prog_1_, true);
    scc4474->add_rule(new parallel_acopy(rel_Prog_1_, rel_Prog_1_1, DELTA, {1, 0}));

    RAM* scc4475 = new RAM(false, 7);
    scc4475->add_relation(rel_App_4_1, true);
    scc4475->add_relation(rel_App_4_1_2_3_4, true);
    scc4475->add_rule(new parallel_acopy(rel_App_4_1, rel_App_4_1_2_3_4, DELTA, {0, 4, 1, 2, 3}));

    RAM* scc4476 = new RAM(false, 2);
    scc4476->add_relation(rel_Lam_4_, true);
    scc4476->add_relation(rel_Lam_4_1_2_3_4, true);
    scc4476->add_rule(new parallel_acopy(rel_Lam_4_, rel_Lam_4_1_2_3_4, DELTA, {4, 0, 1, 2, 3}));

    RAM* scc4477 = new RAM(false, 6);
    scc4477->add_relation(rel_Prog_1_1, false);
    scc4477->add_relation(rel_Time_6_1_2_3_4_5_6, true);
    scc4477->add_rule(new parallel_copy(rel_Time_6_1_2_3_4_5_6, rel_Prog_1_1, FULL, {0, 0, 0, 0, 0, 0}));

    RAM* scc4478 = new RAM(true, 10);
    scc4478->add_relation(rel_AEval_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14, true);
    scc4478->add_relation(rel_Store_14_1, true);
    scc4478->add_relation(rel_App_4_1, false);
    scc4478->add_relation(rel_INT_B_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16, true);
    scc4478->add_relation(rel_INT_B_16_1, true);
    scc4478->add_relation(rel_Var_2_1, false);
    scc4478->add_relation(rel_INT0_17_11, true);
    scc4478->add_relation(rel_Var_2_2, false);
    scc4478->add_relation(rel_ReachesCfg_7_1_2_3_4_5_6_7, true);
    scc4478->add_relation(rel_Step_21_15, true);
    scc4478->add_relation(rel_Step_21_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21, true);
    scc4478->add_relation(rel_Lam_4_1, false);
    scc4478->add_relation(rel_INT_A_15_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15, true);
    scc4478->add_relation(rel_INT2_27_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27, true);
    scc4478->add_relation(rel_Store_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14, true);
    scc4478->add_relation(rel_ReachesCfg_7_1, true);
    scc4478->add_relation(rel_Time_6_1_2_3_4_5_6, true);
    scc4478->add_relation(rel_AEval_14_7_6_5_4_3_2, true);
    scc4478->add_relation(rel_INT00_10_10_9_8_7_6_5_2, true);
    scc4478->add_relation(rel_INT1_20_10_9_8_7_6_5_4, true);
    scc4478->add_relation(rel_INT_D_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16, true);
    scc4478->add_relation(rel_Time_6_, true);
    scc4478->add_relation(rel_INT_C_17_2_10, true);
    scc4478->add_relation(rel_INT_C_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17, true);
    scc4478->add_relation(rel_INT0_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17, true);
    scc4478->add_relation(rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31, true);
    scc4478->add_relation(rel_INT_D_16_10_2, true);
    scc4478->add_relation(rel_Lam_4_, false);
    scc4478->add_relation(rel_INT2_27_10_9_8_7_6_5_3, true);
    scc4478->add_relation(rel_INT_A_15_15_14_13_12_11_10, true);
    scc4478->add_relation(rel_INT00_10_1_2_3_4_5_6_7_8_9_10, true);
    scc4478->add_relation(rel_INT1_20_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20, true);
    scc4478->add_relation(rel_AEval_14_7_6_5_4_3_2_1, true);
    scc4478->add_rule(new parallel_copy(rel_ReachesCfg_7_1_2_3_4_5_6_7, rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31, DELTA, {14, 0, 1, 2, 3, 4, 5}));
    scc4478->add_rule(new parallel_join(rel_INT_B_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16, rel_AEval_14_7_6_5_4_3_2, FULL, rel_INT_A_15_15_14_13_12_11_10, DELTA, {7, 8, 9, 10, 11, 12, 13, 14, 17, 18, 19, 20, 21, 22, 23, 24}));
    scc4478->add_rule(new parallel_acopy(rel_INT0_17_11, rel_INT0_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17, DELTA, {10, 17, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16}));
    scc4478->add_rule(new parallel_join(rel_INT1_20_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20, rel_Lam_4_1, FULL, rel_INT0_17_11, DELTA, {6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 16, 17, 18, 19, 20, 21, 2, 3, 4}));
    scc4478->add_rule(new parallel_acopy(rel_INT2_27_10_9_8_7_6_5_3, rel_INT2_27_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27, DELTA, {9, 8, 7, 6, 5, 4, 2, 27, 0, 1, 3, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26}));
    scc4478->add_rule(new parallel_copy(rel_Time_6_1_2_3_4_5_6, rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31, DELTA, {0, 1, 2, 3, 4, 5}));
    scc4478->add_rule(new parallel_acopy(rel_INT1_20_10_9_8_7_6_5_4, rel_INT1_20_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20, DELTA, {9, 8, 7, 6, 5, 4, 3, 20, 0, 1, 2, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}));
    scc4478->add_rule(new parallel_join(rel_INT_C_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17, rel_Var_2_1, FULL, rel_INT_B_16_1, DELTA, {0, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18}));
    scc4478->add_rule(new parallel_copy_filter(rel_INT_D_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16, rel_INT_C_17_2_10, DELTA, {3, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17}, [](const u64* const data){ return !(data[0] == data[1]); }));
    scc4478->add_rule(new parallel_acopy(rel_INT00_10_10_9_8_7_6_5_2, rel_INT00_10_1_2_3_4_5_6_7_8_9_10, DELTA, {9, 8, 7, 6, 5, 4, 1, 10, 0, 2, 3}));
    scc4478->add_rule(new parallel_join(rel_INT0_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17, rel_INT00_10_10_9_8_7_6_5_2, FULL, rel_AEval_14_7_6_5_4_3_2_1, DELTA, {8, 6, 9, 10, 5, 4, 3, 2, 1, 0, 12, 13, 14, 15, 16, 17, 18}));
    scc4478->add_rule(new parallel_copy(rel_Step_21_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21, rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31, DELTA, {0, 1, 2, 3, 4, 5, 6, 14, 0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13}));
    scc4478->add_rule(new parallel_acopy(rel_Store_14_1, rel_Store_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14, DELTA, {0, 14, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13}));
    scc4478->add_rule(new parallel_join(rel_AEval_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14, rel_Var_2_2, FULL, rel_Store_14_1, DELTA, {2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}));
    scc4478->add_rule(new parallel_acopy(rel_INT_C_17_2_10, rel_INT_C_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17, DELTA, {1, 9, 17, 0, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16}));
    scc4478->add_rule(new parallel_acopy(rel_INT_B_16_1, rel_INT_B_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16, DELTA, {0, 16, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}));
    scc4478->add_rule(new parallel_join(rel_INT_B_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16, rel_AEval_14_7_6_5_4_3_2, DELTA, rel_INT_A_15_15_14_13_12_11_10, FULL, {7, 8, 9, 10, 11, 12, 13, 14, 17, 18, 19, 20, 21, 22, 23, 24}));
    scc4478->add_rule(new parallel_join(rel_INT_B_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16, rel_AEval_14_7_6_5_4_3_2, DELTA, rel_INT_A_15_15_14_13_12_11_10, DELTA, {7, 8, 9, 10, 11, 12, 13, 14, 17, 18, 19, 20, 21, 22, 23, 24}));
    scc4478->add_rule(new parallel_acopy(rel_Time_6_, rel_Time_6_1_2_3_4_5_6, DELTA, {6, 0, 1, 2, 3, 4, 5}));
    scc4478->add_rule(new parallel_join(rel_INT0_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17, rel_INT00_10_10_9_8_7_6_5_2, DELTA, rel_AEval_14_7_6_5_4_3_2_1, FULL, {8, 6, 9, 10, 5, 4, 3, 2, 1, 0, 12, 13, 14, 15, 16, 17, 18}));
    scc4478->add_rule(new parallel_join(rel_INT0_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17, rel_INT00_10_10_9_8_7_6_5_2, DELTA, rel_AEval_14_7_6_5_4_3_2_1, DELTA, {8, 6, 9, 10, 5, 4, 3, 2, 1, 0, 12, 13, 14, 15, 16, 17, 18}));
    scc4478->add_rule(new parallel_acopy(rel_ReachesCfg_7_1, rel_ReachesCfg_7_1_2_3_4_5_6_7, DELTA, {0, 7, 1, 2, 3, 4, 5, 6}));
    scc4478->add_rule(new parallel_copy(rel_Store_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14, rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31, DELTA, {15, 0, 1, 2, 3, 4, 5, 17, 18, 19, 20, 21, 22, 23}));
    scc4478->add_rule(new parallel_acopy(rel_AEval_14_7_6_5_4_3_2_1, rel_AEval_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14, DELTA, {6, 5, 4, 3, 2, 1, 0, 14, 7, 8, 9, 10, 11, 12, 13}));
    scc4478->add_rule(new parallel_acopy(rel_INT_D_16_10_2, rel_INT_D_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16, DELTA, {9, 1, 16, 0, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15}));
    scc4478->add_rule(new parallel_join(rel_INT2_27_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27, rel_AEval_14_7_6_5_4_3_2_1, DELTA, rel_INT1_20_10_9_8_7_6_5_4, FULL, {16, 17, 18, 6, 5, 4, 3, 2, 1, 0, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 8, 9, 10, 11, 12, 13, 14}));
    scc4478->add_rule(new parallel_join(rel_INT2_27_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27, rel_AEval_14_7_6_5_4_3_2_1, DELTA, rel_INT1_20_10_9_8_7_6_5_4, DELTA, {16, 17, 18, 6, 5, 4, 3, 2, 1, 0, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 8, 9, 10, 11, 12, 13, 14}));
    scc4478->add_rule(new parallel_copy_filter(rel_AEval_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14, rel_INT_D_16_10_2, DELTA, {3, 11, 12, 13, 14, 15, 16, 4, 5, 6, 7, 8, 9, 10}, [](const u64* const data){ return !(data[0] == data[1]); }));
    scc4478->add_rule(new parallel_join(rel_INT_A_15_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15, rel_Lam_4_1, FULL, rel_Step_21_15, DELTA, {0, 2, 3, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25}));
    scc4478->add_rule(new parallel_copy(rel_Store_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14, rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31, DELTA, {16, 0, 1, 2, 3, 4, 5, 24, 25, 26, 27, 28, 29, 30}));
    scc4478->add_rule(new parallel_join(rel_AEval_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14, rel_Lam_4_, FULL, rel_Time_6_, DELTA, {1, 6, 7, 8, 9, 10, 11, 1, 6, 7, 8, 9, 10, 11}));
    scc4478->add_rule(new parallel_join(rel_INT2_27_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27, rel_AEval_14_7_6_5_4_3_2_1, FULL, rel_INT1_20_10_9_8_7_6_5_4, DELTA, {16, 17, 18, 6, 5, 4, 3, 2, 1, 0, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 8, 9, 10, 11, 12, 13, 14}));
    scc4478->add_rule(new parallel_acopy(rel_AEval_14_7_6_5_4_3_2, rel_AEval_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14, DELTA, {6, 5, 4, 3, 2, 1, 14, 0, 7, 8, 9, 10, 11, 12, 13}));
    scc4478->add_rule(new parallel_acopy(rel_Step_21_15, rel_Step_21_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21, DELTA, {14, 21, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20}));
    scc4478->add_rule(new parallel_join(rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31, rel_AEval_14_7_6_5_4_3_2_1, DELTA, rel_INT2_27_10_9_8_7_6_5_3, FULL, {16, 5, 4, 3, 2, 1, 0, 19, 20, 21, 22, 23, 24, 25, 28, 26, 27, 8, 9, 10, 11, 12, 13, 14, 29, 30, 31, 32, 33, 34, 35}));
    scc4478->add_rule(new parallel_join(rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31, rel_AEval_14_7_6_5_4_3_2_1, FULL, rel_INT2_27_10_9_8_7_6_5_3, DELTA, {16, 5, 4, 3, 2, 1, 0, 19, 20, 21, 22, 23, 24, 25, 28, 26, 27, 8, 9, 10, 11, 12, 13, 14, 29, 30, 31, 32, 33, 34, 35}));
    scc4478->add_rule(new parallel_join(rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31, rel_AEval_14_7_6_5_4_3_2_1, DELTA, rel_INT2_27_10_9_8_7_6_5_3, DELTA, {16, 5, 4, 3, 2, 1, 0, 19, 20, 21, 22, 23, 24, 25, 28, 26, 27, 8, 9, 10, 11, 12, 13, 14, 29, 30, 31, 32, 33, 34, 35}));
    scc4478->add_rule(new parallel_join(rel_INT00_10_1_2_3_4_5_6_7_8_9_10, rel_App_4_1, FULL, rel_ReachesCfg_7_1, DELTA, {0, 2, 3, 4, 6, 7, 8, 9, 10, 11}));
    scc4478->add_rule(new parallel_acopy(rel_INT_A_15_15_14_13_12_11_10, rel_INT_A_15_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15, DELTA, {14, 13, 12, 11, 10, 9, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8}));

    RAM* scc4479 = new RAM(false, 4);
    scc4479->add_relation(rel_ReachesCfg_7_1_2_3_4_5_6_7, true);
    scc4479->add_relation(rel_Prog_1_1, false);
    scc4479->add_rule(new parallel_copy(rel_ReachesCfg_7_1_2_3_4_5_6_7, rel_Prog_1_1, FULL, {0, 0, 0, 0, 0, 0, 0}));

    RAM* scc4480 = new RAM(false, 8);
    scc4480->add_relation(rel_Lam_4_1, true);
    scc4480->add_relation(rel_Lam_4_1_2_3_4, true);
    scc4480->add_rule(new parallel_acopy(rel_Lam_4_1, rel_Lam_4_1_2_3_4, DELTA, {0, 4, 1, 2, 3}));

    LIE* lie = new LIE();
    lie->add_relation(rel_Lam_4_1_2_3_4);
    lie->add_relation(rel_AEval_14_7_6_5_4_3_2_1);
    lie->add_relation(rel_App_4_1_2_3_4);
    lie->add_relation(rel_INT1_20_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20);
    lie->add_relation(rel_INT00_10_1_2_3_4_5_6_7_8_9_10);
    lie->add_relation(rel_INT_A_15_15_14_13_12_11_10);
    lie->add_relation(rel_INT2_27_10_9_8_7_6_5_3);
    lie->add_relation(rel_Lam_4_);
    lie->add_relation(rel_INT_D_16_10_2);
    lie->add_relation(rel_Prog_1_);
    lie->add_relation(rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31);
    lie->add_relation(rel_INT0_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17);
    lie->add_relation(rel_INT_C_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17);
    lie->add_relation(rel_INT_C_17_2_10);
    lie->add_relation(rel_Time_6_);
    lie->add_relation(rel_INT_D_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16);
    lie->add_relation(rel_INT1_20_10_9_8_7_6_5_4);
    lie->add_relation(rel_INT00_10_10_9_8_7_6_5_2);
    lie->add_relation(rel_AEval_14_7_6_5_4_3_2);
    lie->add_relation(rel_ReachesClo_7_1_2_3_4_5_6_7);
    lie->add_relation(rel_Time_6_1_2_3_4_5_6);
    lie->add_relation(rel_ReachesCfg_7_1);
    lie->add_relation(rel_Store_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14);
    lie->add_relation(rel_INT2_27_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27);
    lie->add_relation(rel_INT_A_15_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15);
    lie->add_relation(rel_Lam_4_1);
    lie->add_relation(rel_Step_21_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21);
    lie->add_relation(rel_Step_21_15);
    lie->add_relation(rel_Prog_1_1);
    lie->add_relation(rel_ReachesCfg_7_1_2_3_4_5_6_7);
    lie->add_relation(rel_Var_2_2);
    lie->add_relation(rel_INT0_17_11);
    lie->add_relation(rel_Var_2_1);
    lie->add_relation(rel_Var_2_1_2);
    lie->add_relation(rel_INT_B_16_1);
    lie->add_relation(rel_INT_B_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16);
    lie->add_relation(rel_App_4_1);
    lie->add_relation(rel_Store_14_1);
    lie->add_relation(rel_AEval_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14);
    lie->add_scc(scc4471);
    lie->add_scc(scc4472);
    lie->add_scc(scc4473);
    lie->add_scc(scc4474);
    lie->add_scc(scc4475);
    lie->add_scc(scc4476);
    lie->add_scc(scc4477);
    lie->add_scc(scc4478);
    lie->add_scc(scc4479);
    lie->add_scc(scc4480);
    lie->add_scc_dependance(scc4471, scc4478);
    lie->add_scc_dependance(scc4473, scc4478);
    lie->add_scc_dependance(scc4475, scc4478);
    lie->add_scc_dependance(scc4476, scc4478);
    lie->add_scc_dependance(scc4477, scc4478);
    lie->add_scc_dependance(scc4478, scc4472);
    lie->add_scc_dependance(scc4479, scc4478);
    lie->add_scc_dependance(scc4480, scc4478);



    lie->set_comm(mcomm);
    lie->set_batch_size(1);
    lie->execute();

    delete lie;

    mcomm.destroy();
    return 0;
}
