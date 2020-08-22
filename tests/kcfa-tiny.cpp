#include "../src/parallel_RA_inc.h"


int main(int argc, char **argv)
{
    mpi_comm mcomm;
    mcomm.create(argc, argv);


    relation* rel_Free0_4_2_1 = new relation(2, false, 4, 258, "rel_Free0_4_2_1", "../data/g4470/Free0_4_40", FULL);
    relation* rel_Lam_4_1_2_3_4 = new relation(4, true, 4, 262, "rel_Lam_4_1_2_3_4", "../data/g4470/Lam_4_39", FULL);
    relation* rel_AEval_14_7_6_5_4_3_2_1 = new relation(7, false, 14, 275, "rel_AEval_14_7_6_5_4_3_2_1", "../data/g4470/AEval_14_38", FULL);
    relation* rel_App_4_1_2_3_4 = new relation(4, true, 4, 273, "rel_App_4_1_2_3_4", "../data/g4470/App_4_37", FULL);
    relation* rel_INT1_20_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20 = new relation(20, true, 20, 261, "rel_INT1_20_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20", "../data/g4470/INT1_20_36", FULL);
    relation* rel_Free_2_2 = new relation(1, false, 2, 264, "rel_Free_2_2", "../data/g4470/Free_2_35", FULL);
    relation* rel_INT00_10_1_2_3_4_5_6_7_8_9_10 = new relation(10, true, 10, 266, "rel_INT00_10_1_2_3_4_5_6_7_8_9_10", "../data/g4470/INT00_10_34", FULL);
    relation* rel_Store_14_7_6_5_4_3_2_1 = new relation(7, false, 14, 272, "rel_Store_14_7_6_5_4_3_2_1", "../data/g4470/Store_14_33", FULL);
    relation* rel_FrProp_13_12_11_10_9_8_7_13 = new relation(7, false, 13, 268, "rel_FrProp_13_12_11_10_9_8_7_13", "../data/g4470/FrProp_13_32", FULL);
    relation* rel_INT2_27_10_9_8_7_6_5_3 = new relation(7, false, 27, 270, "rel_INT2_27_10_9_8_7_6_5_3", "../data/g4470/INT2_27_31", FULL);
    relation* rel_Prog_1_ = new relation(0, false, 1, 259, "rel_Prog_1_", "../data/g4470/Prog_1_30", FULL);
    relation* rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31 = new relation(31, true, 31, 269, "rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31", "../data/g4470/APP_31_29", FULL);
    relation* rel_Free1_3_1_2_3 = new relation(3, true, 3, 260, "rel_Free1_3_1_2_3", "../data/g4470/Free1_3_28", FULL);
    relation* rel_INT0_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17 = new relation(17, true, 17, 256, "rel_INT0_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17", "../data/g4470/INT0_17_27", FULL);
    relation* rel_Free1_3_2_1 = new relation(2, false, 3, 260, "rel_Free1_3_2_1", "../data/g4470/Free1_3_26", FULL);
    relation* rel_INT1_20_10_9_8_7_6_5_4 = new relation(7, false, 20, 261, "rel_INT1_20_10_9_8_7_6_5_4", "../data/g4470/INT1_20_25", FULL);
    relation* rel_INT00_10_10_9_8_7_6_5_2 = new relation(7, false, 10, 266, "rel_INT00_10_10_9_8_7_6_5_2", "../data/g4470/INT00_10_24", FULL);
    relation* rel_ReachesClo_7_1_2_3_4_5_6_7 = new relation(7, true, 7, 265, "rel_ReachesClo_7_1_2_3_4_5_6_7", "../data/g4470/ReachesClo_7_23", FULL);
    relation* rel_App_4_2 = new relation(1, false, 4, 273, "rel_App_4_2", "../data/g4470/App_4_22", FULL);
    relation* rel_Time_6_1_2_3_4_5_6 = new relation(6, true, 6, 274, "rel_Time_6_1_2_3_4_5_6", "../data/g4470/Time_6_21", FULL);
    relation* rel_ReachesCfg_7_1 = new relation(1, false, 7, 257, "rel_ReachesCfg_7_1", "../data/g4470/ReachesCfg_7_20", FULL);
    relation* rel_App_4_3 = new relation(1, false, 4, 273, "rel_App_4_3", "../data/g4470/App_4_19", FULL);
    relation* rel_AE0_8_8_7_6_5_4_3_2 = new relation(7, false, 8, 263, "rel_AE0_8_8_7_6_5_4_3_2", "../data/g4470/AE0_8_18", FULL);
    relation* rel_Store_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14 = new relation(14, true, 14, 272, "rel_Store_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14", "../data/g4470/Store_14_17", FULL);
    relation* rel_INT2_27_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27 = new relation(27, true, 27, 270, "rel_INT2_27_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27", "../data/g4470/INT2_27_16", FULL);
    relation* rel_Lam_4_4 = new relation(1, false, 4, 262, "rel_Lam_4_4", "../data/g4470/Lam_4_15", FULL);
    relation* rel_Lam_4_1 = new relation(1, false, 4, 262, "rel_Lam_4_1", "../data/g4470/Lam_4_14", FULL);
    relation* rel_Step_21_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21 = new relation(21, true, 21, 276, "rel_Step_21_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21", "../data/g4470/Step_21_13", FULL);
    relation* rel_Step_21_15 = new relation(1, false, 21, 276, "rel_Step_21_15", "../data/g4470/Step_21_12", FULL);
    relation* rel_App_4_4 = new relation(1, false, 4, 273, "rel_App_4_4", "../data/g4470/App_4_11", FULL);
    relation* rel_Prog_1_1 = new relation(1, true, 1, 259, "rel_Prog_1_1", "../data/g4470/Prog_1_10", FULL);
    relation* rel_ReachesCfg_7_1_2_3_4_5_6_7 = new relation(7, true, 7, 257, "rel_ReachesCfg_7_1_2_3_4_5_6_7", "../data/g4470/ReachesCfg_7_9", FULL);
    relation* rel_INT0_17_11 = new relation(1, false, 17, 256, "rel_INT0_17_11", "../data/g4470/INT0_17_8", FULL);
    relation* rel_Free0_4_1_2_3_4 = new relation(4, true, 4, 258, "rel_Free0_4_1_2_3_4", "../data/g4470/Free0_4_7", FULL);
    relation* rel_FrProp_13_1_2_3_4_5_6_7_8_9_10_11_12_13 = new relation(13, true, 13, 268, "rel_FrProp_13_1_2_3_4_5_6_7_8_9_10_11_12_13", "../data/g4470/FrProp_13_6", FULL);
    relation* rel_Var_2_1 = new relation(1, false, 2, 267, "rel_Var_2_1", "../data/g4470/Var_2_5", FULL);
    relation* rel_Free_2_1_2 = new relation(2, true, 2, 264, "rel_Free_2_1_2", "../data/g4470/Free_2_4", FULL);
    relation* rel_Var_2_1_2 = new relation(2, true, 2, 267, "rel_Var_2_1_2", "../data/g4470/Var_2_3", FULL);
    relation* rel_AE0_8_1_2_3_4_5_6_7_8 = new relation(8, true, 8, 263, "rel_AE0_8_1_2_3_4_5_6_7_8", "../data/g4470/AE0_8_2", FULL);
    relation* rel_App_4_1 = new relation(1, false, 4, 273, "rel_App_4_1", "../data/g4470/App_4_1", FULL);
    relation* rel_AEval_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14 = new relation(14, true, 14, 275, "rel_AEval_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14", "../data/g4470/AEval_14_0", FULL);

    RAM* scc4471 = new RAM(false, 1);
    scc4471->add_relation(rel_Time_6_1_2_3_4_5_6, true);
    scc4471->add_relation(rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31, false);
    scc4471->add_rule(new parallel_copy(rel_Time_6_1_2_3_4_5_6, rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31, FULL, {0, 1, 2, 3, 4, 5}));

    RAM* scc4472 = new RAM(false, 5);
    scc4472->add_relation(rel_App_4_4, true);
    scc4472->add_relation(rel_App_4_1_2_3_4, true);
    scc4472->add_rule(new parallel_acopy(rel_App_4_4, rel_App_4_1_2_3_4, DELTA, {3, 4, 0, 1, 2}));

    RAM* scc4473 = new RAM(false, 9);
    scc4473->add_relation(rel_Var_2_1_2, false);
    scc4473->add_relation(rel_Free_2_1_2, true);
    scc4473->add_rule(new parallel_copy(rel_Free_2_1_2, rel_Var_2_1_2, FULL, {1, 0}));

    RAM* scc4474 = new RAM(false, 13);
    scc4474->add_relation(rel_Lam_4_1, true);
    scc4474->add_relation(rel_Lam_4_1_2_3_4, true);
    scc4474->add_rule(new parallel_acopy(rel_Lam_4_1, rel_Lam_4_1_2_3_4, DELTA, {0, 4, 1, 2, 3}));

    RAM* scc4475 = new RAM(false, 3);
    scc4475->add_relation(rel_ReachesClo_7_1_2_3_4_5_6_7, true);
    scc4475->add_relation(rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31, false);
    scc4475->add_rule(new parallel_copy(rel_ReachesClo_7_1_2_3_4_5_6_7, rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31, FULL, {7, 8, 9, 10, 11, 12, 13}));

    RAM* scc4476 = new RAM(true, 7);
    scc4476->add_relation(rel_AEval_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14, true);
    scc4476->add_relation(rel_App_4_1, false);
    scc4476->add_relation(rel_AE0_8_1_2_3_4_5_6_7_8, true);
    scc4476->add_relation(rel_Var_2_1, false);
    scc4476->add_relation(rel_FrProp_13_1_2_3_4_5_6_7_8_9_10_11_12_13, true);
    scc4476->add_relation(rel_INT0_17_11, true);
    scc4476->add_relation(rel_ReachesCfg_7_1_2_3_4_5_6_7, true);
    scc4476->add_relation(rel_Step_21_15, true);
    scc4476->add_relation(rel_Step_21_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21, true);
    scc4476->add_relation(rel_Lam_4_1, false);
    scc4476->add_relation(rel_INT2_27_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27, true);
    scc4476->add_relation(rel_Store_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14, true);
    scc4476->add_relation(rel_AE0_8_8_7_6_5_4_3_2, true);
    scc4476->add_relation(rel_ReachesCfg_7_1, true);
    scc4476->add_relation(rel_INT00_10_10_9_8_7_6_5_2, true);
    scc4476->add_relation(rel_INT1_20_10_9_8_7_6_5_4, true);
    scc4476->add_relation(rel_INT0_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17, true);
    scc4476->add_relation(rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31, true);
    scc4476->add_relation(rel_INT2_27_10_9_8_7_6_5_3, true);
    scc4476->add_relation(rel_FrProp_13_12_11_10_9_8_7_13, true);
    scc4476->add_relation(rel_Store_14_7_6_5_4_3_2_1, true);
    scc4476->add_relation(rel_INT00_10_1_2_3_4_5_6_7_8_9_10, true);
    scc4476->add_relation(rel_Free_2_2, false);
    scc4476->add_relation(rel_INT1_20_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20, true);
    scc4476->add_relation(rel_AEval_14_7_6_5_4_3_2_1, true);
    scc4476->add_rule(new parallel_copy(rel_Step_21_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21, rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31, DELTA, {0, 1, 2, 3, 4, 5, 6, 14, 0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13}));
    scc4476->add_rule(new parallel_join(rel_Store_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14, rel_FrProp_13_12_11_10_9_8_7_13, DELTA, rel_Store_14_7_6_5_4_3_2_1, FULL, {6, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 21}));
    scc4476->add_rule(new parallel_join(rel_Store_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14, rel_FrProp_13_12_11_10_9_8_7_13, DELTA, rel_Store_14_7_6_5_4_3_2_1, DELTA, {6, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 21}));
    scc4476->add_rule(new parallel_acopy(rel_INT0_17_11, rel_INT0_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17, DELTA, {10, 17, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16}));
    scc4476->add_rule(new parallel_join(rel_INT1_20_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20, rel_Lam_4_1, FULL, rel_INT0_17_11, DELTA, {6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 16, 17, 18, 19, 20, 21, 2, 3, 4}));
    scc4476->add_rule(new parallel_copy(rel_Store_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14, rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31, DELTA, {15, 0, 1, 2, 3, 4, 5, 17, 18, 19, 20, 21, 22, 23}));
    scc4476->add_rule(new parallel_join(rel_AEval_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14, rel_AE0_8_8_7_6_5_4_3_2, DELTA, rel_Store_14_7_6_5_4_3_2_1, FULL, {8, 5, 4, 3, 2, 1, 0, 10, 11, 12, 13, 14, 15, 16}));
    scc4476->add_rule(new parallel_acopy(rel_INT2_27_10_9_8_7_6_5_3, rel_INT2_27_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27, DELTA, {9, 8, 7, 6, 5, 4, 2, 27, 0, 1, 3, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26}));
    scc4476->add_rule(new parallel_acopy(rel_INT1_20_10_9_8_7_6_5_4, rel_INT1_20_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20, DELTA, {9, 8, 7, 6, 5, 4, 3, 20, 0, 1, 2, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}));
    scc4476->add_rule(new parallel_copy(rel_ReachesCfg_7_1_2_3_4_5_6_7, rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31, DELTA, {14, 0, 1, 2, 3, 4, 5}));
    scc4476->add_rule(new parallel_join(rel_AEval_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14, rel_AE0_8_8_7_6_5_4_3_2, FULL, rel_Store_14_7_6_5_4_3_2_1, DELTA, {8, 5, 4, 3, 2, 1, 0, 10, 11, 12, 13, 14, 15, 16}));
    scc4476->add_rule(new parallel_join(rel_AEval_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14, rel_Lam_4_1, FULL, rel_ReachesCfg_7_1, DELTA, {0, 6, 7, 8, 9, 10, 11, 0, 6, 7, 8, 9, 10, 11}));
    scc4476->add_rule(new parallel_acopy(rel_INT00_10_10_9_8_7_6_5_2, rel_INT00_10_1_2_3_4_5_6_7_8_9_10, DELTA, {9, 8, 7, 6, 5, 4, 1, 10, 0, 2, 3}));
    scc4476->add_rule(new parallel_join(rel_AE0_8_1_2_3_4_5_6_7_8, rel_Var_2_1, FULL, rel_ReachesCfg_7_1, DELTA, {0, 2, 4, 5, 6, 7, 8, 9}));
    scc4476->add_rule(new parallel_join(rel_INT2_27_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27, rel_AEval_14_7_6_5_4_3_2_1, FULL, rel_INT1_20_10_9_8_7_6_5_4, DELTA, {16, 17, 18, 6, 5, 4, 3, 2, 1, 0, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 8, 9, 10, 11, 12, 13, 14}));
    scc4476->add_rule(new parallel_join(rel_INT0_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17, rel_INT00_10_10_9_8_7_6_5_2, FULL, rel_AEval_14_7_6_5_4_3_2_1, DELTA, {8, 6, 9, 10, 5, 4, 3, 2, 1, 0, 12, 13, 14, 15, 16, 17, 18}));
    scc4476->add_rule(new parallel_join(rel_ReachesCfg_7_1_2_3_4_5_6_7, rel_App_4_1, FULL, rel_ReachesCfg_7_1, DELTA, {4, 6, 7, 8, 9, 10, 11}));
    scc4476->add_rule(new parallel_join(rel_AEval_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14, rel_AE0_8_8_7_6_5_4_3_2, DELTA, rel_Store_14_7_6_5_4_3_2_1, DELTA, {8, 5, 4, 3, 2, 1, 0, 10, 11, 12, 13, 14, 15, 16}));
    scc4476->add_rule(new parallel_acopy(rel_FrProp_13_12_11_10_9_8_7_13, rel_FrProp_13_1_2_3_4_5_6_7_8_9_10_11_12_13, DELTA, {11, 10, 9, 8, 7, 6, 12, 13, 0, 1, 2, 3, 4, 5}));
    scc4476->add_rule(new parallel_copy(rel_Store_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14, rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31, DELTA, {16, 0, 1, 2, 3, 4, 5, 24, 25, 26, 27, 28, 29, 30}));
    scc4476->add_rule(new parallel_join(rel_INT2_27_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27, rel_AEval_14_7_6_5_4_3_2_1, DELTA, rel_INT1_20_10_9_8_7_6_5_4, FULL, {16, 17, 18, 6, 5, 4, 3, 2, 1, 0, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 8, 9, 10, 11, 12, 13, 14}));
    scc4476->add_rule(new parallel_join(rel_INT2_27_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27, rel_AEval_14_7_6_5_4_3_2_1, DELTA, rel_INT1_20_10_9_8_7_6_5_4, DELTA, {16, 17, 18, 6, 5, 4, 3, 2, 1, 0, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 8, 9, 10, 11, 12, 13, 14}));
    scc4476->add_rule(new parallel_acopy(rel_AE0_8_8_7_6_5_4_3_2, rel_AE0_8_1_2_3_4_5_6_7_8, DELTA, {7, 6, 5, 4, 3, 2, 1, 8, 0}));
    scc4476->add_rule(new parallel_join(rel_Store_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14, rel_FrProp_13_12_11_10_9_8_7_13, FULL, rel_Store_14_7_6_5_4_3_2_1, DELTA, {6, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 21}));
    scc4476->add_rule(new parallel_acopy(rel_ReachesCfg_7_1, rel_ReachesCfg_7_1_2_3_4_5_6_7, DELTA, {0, 7, 1, 2, 3, 4, 5, 6}));
    scc4476->add_rule(new parallel_join(rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31, rel_AEval_14_7_6_5_4_3_2_1, DELTA, rel_INT2_27_10_9_8_7_6_5_3, FULL, {16, 5, 4, 3, 2, 1, 0, 19, 20, 21, 22, 23, 24, 25, 28, 26, 27, 8, 9, 10, 11, 12, 13, 14, 29, 30, 31, 32, 33, 34, 35}));
    scc4476->add_rule(new parallel_join(rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31, rel_AEval_14_7_6_5_4_3_2_1, FULL, rel_INT2_27_10_9_8_7_6_5_3, DELTA, {16, 5, 4, 3, 2, 1, 0, 19, 20, 21, 22, 23, 24, 25, 28, 26, 27, 8, 9, 10, 11, 12, 13, 14, 29, 30, 31, 32, 33, 34, 35}));
    scc4476->add_rule(new parallel_join(rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31, rel_AEval_14_7_6_5_4_3_2_1, DELTA, rel_INT2_27_10_9_8_7_6_5_3, DELTA, {16, 5, 4, 3, 2, 1, 0, 19, 20, 21, 22, 23, 24, 25, 28, 26, 27, 8, 9, 10, 11, 12, 13, 14, 29, 30, 31, 32, 33, 34, 35}));
    scc4476->add_rule(new parallel_join(rel_INT0_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17, rel_INT00_10_10_9_8_7_6_5_2, DELTA, rel_AEval_14_7_6_5_4_3_2_1, FULL, {8, 6, 9, 10, 5, 4, 3, 2, 1, 0, 12, 13, 14, 15, 16, 17, 18}));
    scc4476->add_rule(new parallel_join(rel_INT0_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17, rel_INT00_10_10_9_8_7_6_5_2, DELTA, rel_AEval_14_7_6_5_4_3_2_1, DELTA, {8, 6, 9, 10, 5, 4, 3, 2, 1, 0, 12, 13, 14, 15, 16, 17, 18}));
    scc4476->add_rule(new parallel_acopy(rel_AEval_14_7_6_5_4_3_2_1, rel_AEval_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14, DELTA, {6, 5, 4, 3, 2, 1, 0, 14, 7, 8, 9, 10, 11, 12, 13}));
    scc4476->add_rule(new parallel_join(rel_ReachesCfg_7_1_2_3_4_5_6_7, rel_App_4_1, FULL, rel_ReachesCfg_7_1, DELTA, {3, 6, 7, 8, 9, 10, 11}));
    scc4476->add_rule(new parallel_acopy(rel_Store_14_7_6_5_4_3_2_1, rel_Store_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14, DELTA, {6, 5, 4, 3, 2, 1, 0, 14, 7, 8, 9, 10, 11, 12, 13}));
    scc4476->add_rule(new parallel_join(rel_INT00_10_1_2_3_4_5_6_7_8_9_10, rel_App_4_1, FULL, rel_ReachesCfg_7_1, DELTA, {0, 2, 3, 4, 6, 7, 8, 9, 10, 11}));
    scc4476->add_rule(new parallel_join(rel_FrProp_13_1_2_3_4_5_6_7_8_9_10_11_12_13, rel_Free_2_2, FULL, rel_Step_21_15, DELTA, {12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 2}));
    scc4476->add_rule(new parallel_acopy(rel_Step_21_15, rel_Step_21_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21, DELTA, {14, 21, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20}));
    scc4476->add_rule(new parallel_join(rel_ReachesCfg_7_1_2_3_4_5_6_7, rel_App_4_1, FULL, rel_ReachesCfg_7_1, DELTA, {2, 6, 7, 8, 9, 10, 11}));

    RAM* scc4477 = new RAM(false, 11);
    scc4477->add_relation(rel_Prog_1_1, false);
    scc4477->add_relation(rel_Time_6_1_2_3_4_5_6, true);
    scc4477->add_rule(new parallel_copy(rel_Time_6_1_2_3_4_5_6, rel_Prog_1_1, FULL, {0, 0, 0, 0, 0, 0}));

    RAM* scc4478 = new RAM(false, 15);
    scc4478->add_relation(rel_App_4_2, true);
    scc4478->add_relation(rel_App_4_1_2_3_4, true);
    scc4478->add_rule(new parallel_acopy(rel_App_4_2, rel_App_4_1_2_3_4, DELTA, {1, 4, 0, 2, 3}));

    RAM* scc4479 = new RAM(false, 2);
    scc4479->add_relation(rel_Lam_4_4, true);
    scc4479->add_relation(rel_Lam_4_1_2_3_4, true);
    scc4479->add_rule(new parallel_acopy(rel_Lam_4_4, rel_Lam_4_1_2_3_4, DELTA, {3, 4, 0, 1, 2}));

    RAM* scc4480 = new RAM(false, 6);
    scc4480->add_relation(rel_Prog_1_1, true);
    scc4480->add_relation(rel_Prog_1_, true);
    scc4480->add_rule(new parallel_acopy(rel_Prog_1_, rel_Prog_1_1, DELTA, {1, 0}));

    RAM* scc4481 = new RAM(true, 10);
    scc4481->add_relation(rel_Free_2_1_2, true);
    scc4481->add_relation(rel_Free0_4_1_2_3_4, true);
    scc4481->add_relation(rel_App_4_4, false);
    scc4481->add_relation(rel_Lam_4_4, false);
    scc4481->add_relation(rel_App_4_3, false);
    scc4481->add_relation(rel_App_4_2, false);
    scc4481->add_relation(rel_Free1_3_2_1, true);
    scc4481->add_relation(rel_Free1_3_1_2_3, true);
    scc4481->add_relation(rel_Free_2_2, true);
    scc4481->add_relation(rel_Free0_4_2_1, true);
    scc4481->add_rule(new parallel_copy_filter(rel_Free1_3_1_2_3, rel_Free0_4_2_1, DELTA, {1, 3, 4}, [](const u64* const data){ return !(data[0] == data[1]); }));
    scc4481->add_rule(new parallel_acopy(rel_Free0_4_2_1, rel_Free0_4_1_2_3_4, DELTA, {1, 0, 4, 2, 3}));
    scc4481->add_rule(new parallel_join(rel_Free_2_1_2, rel_Free_2_2, DELTA, rel_App_4_4, FULL, {2, 4}));
    scc4481->add_rule(new parallel_join(rel_Free0_4_1_2_3_4, rel_Free_2_2, DELTA, rel_Lam_4_4, FULL, {2, 5, 6, 4}));
    scc4481->add_rule(new parallel_copy_filter(rel_Free_2_1_2, rel_Free1_3_2_1, DELTA, {1, 3}, [](const u64* const data){ return !(data[0] == data[1]); }));
    scc4481->add_rule(new parallel_join(rel_Free_2_1_2, rel_Free_2_2, DELTA, rel_App_4_3, FULL, {2, 4}));
    scc4481->add_rule(new parallel_acopy(rel_Free1_3_2_1, rel_Free1_3_1_2_3, DELTA, {1, 0, 3, 2}));
    scc4481->add_rule(new parallel_acopy(rel_Free_2_2, rel_Free_2_1_2, DELTA, {1, 2, 0}));
    scc4481->add_rule(new parallel_join(rel_Free_2_1_2, rel_Free_2_2, DELTA, rel_App_4_2, FULL, {2, 4}));

    RAM* scc4482 = new RAM(false, 14);
    scc4482->add_relation(rel_App_4_3, true);
    scc4482->add_relation(rel_App_4_1_2_3_4, true);
    scc4482->add_rule(new parallel_acopy(rel_App_4_3, rel_App_4_1_2_3_4, DELTA, {2, 4, 0, 1, 3}));

    RAM* scc4483 = new RAM(false, 4);
    scc4483->add_relation(rel_Var_2_1_2, true);
    scc4483->add_relation(rel_Var_2_1, true);
    scc4483->add_rule(new parallel_acopy(rel_Var_2_1, rel_Var_2_1_2, DELTA, {0, 2, 1}));

    RAM* scc4484 = new RAM(false, 8);
    scc4484->add_relation(rel_ReachesCfg_7_1_2_3_4_5_6_7, true);
    scc4484->add_relation(rel_Prog_1_1, false);
    scc4484->add_rule(new parallel_copy(rel_ReachesCfg_7_1_2_3_4_5_6_7, rel_Prog_1_1, FULL, {0, 0, 0, 0, 0, 0, 0}));

    RAM* scc4485 = new RAM(false, 12);
    scc4485->add_relation(rel_App_4_1, true);
    scc4485->add_relation(rel_App_4_1_2_3_4, true);
    scc4485->add_rule(new parallel_acopy(rel_App_4_1, rel_App_4_1_2_3_4, DELTA, {0, 4, 1, 2, 3}));

    LIE* lie = new LIE();
    lie->add_relation(rel_Free0_4_2_1);
    lie->add_relation(rel_Lam_4_1_2_3_4);
    lie->add_relation(rel_AEval_14_7_6_5_4_3_2_1);
    lie->add_relation(rel_App_4_1_2_3_4);
    lie->add_relation(rel_INT1_20_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20);
    lie->add_relation(rel_Free_2_2);
    lie->add_relation(rel_INT00_10_1_2_3_4_5_6_7_8_9_10);
    lie->add_relation(rel_Store_14_7_6_5_4_3_2_1);
    lie->add_relation(rel_FrProp_13_12_11_10_9_8_7_13);
    lie->add_relation(rel_INT2_27_10_9_8_7_6_5_3);
    lie->add_relation(rel_Prog_1_);
    lie->add_relation(rel_APP_31_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31);
    lie->add_relation(rel_Free1_3_1_2_3);
    lie->add_relation(rel_INT0_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17);
    lie->add_relation(rel_Free1_3_2_1);
    lie->add_relation(rel_INT1_20_10_9_8_7_6_5_4);
    lie->add_relation(rel_INT00_10_10_9_8_7_6_5_2);
    lie->add_relation(rel_ReachesClo_7_1_2_3_4_5_6_7);
    lie->add_relation(rel_App_4_2);
    lie->add_relation(rel_Time_6_1_2_3_4_5_6);
    lie->add_relation(rel_ReachesCfg_7_1);
    lie->add_relation(rel_App_4_3);
    lie->add_relation(rel_AE0_8_8_7_6_5_4_3_2);
    lie->add_relation(rel_Store_14_1_2_3_4_5_6_7_8_9_10_11_12_13_14);
    lie->add_relation(rel_INT2_27_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27);
    lie->add_relation(rel_Lam_4_4);
    lie->add_relation(rel_Lam_4_1);
    lie->add_relation(rel_Step_21_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21);
    lie->add_relation(rel_Step_21_15);
    lie->add_relation(rel_App_4_4);
    lie->add_relation(rel_Prog_1_1);
    lie->add_relation(rel_ReachesCfg_7_1_2_3_4_5_6_7);
    lie->add_relation(rel_INT0_17_11);
    lie->add_relation(rel_Free0_4_1_2_3_4);
    lie->add_relation(rel_FrProp_13_1_2_3_4_5_6_7_8_9_10_11_12_13);
    lie->add_relation(rel_Var_2_1);
    lie->add_relation(rel_Free_2_1_2);
    lie->add_relation(rel_Var_2_1_2);
    lie->add_relation(rel_AE0_8_1_2_3_4_5_6_7_8);
    lie->add_relation(rel_App_4_1);
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
    lie->add_scc(scc4481);
    lie->add_scc(scc4482);
    lie->add_scc(scc4483);
    lie->add_scc(scc4484);
    lie->add_scc(scc4485);
    lie->add_scc_dependance(scc4472, scc4481);
    lie->add_scc_dependance(scc4473, scc4481);
    lie->add_scc_dependance(scc4474, scc4476);
    lie->add_scc_dependance(scc4476, scc4475);
    lie->add_scc_dependance(scc4476, scc4471);
    lie->add_scc_dependance(scc4478, scc4481);
    lie->add_scc_dependance(scc4479, scc4481);
    lie->add_scc_dependance(scc4481, scc4476);
    lie->add_scc_dependance(scc4482, scc4481);
    lie->add_scc_dependance(scc4483, scc4476);
    lie->add_scc_dependance(scc4484, scc4476);
    lie->add_scc_dependance(scc4485, scc4476);




    lie->set_comm(mcomm);
    lie->set_batch_size(1);
    lie->execute();

    delete lie;

    mcomm.destroy();
    return 0;
}
