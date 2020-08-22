#include "../src/parallel_RA_inc.h"


int main(int argc, char **argv)
{
    mpi_comm mcomm;
    mcomm.create(argc, argv);


    relation* rel_Free0_4_2_1 = new relation(2, false, 4, 258, "rel_Free0_4_2_1", "../data/g4470/Free0_4_40", FULL);
    relation* rel_Lam_4_1_2_3_4 = new relation(4, true, 4, 264, "rel_Lam_4_1_2_3_4", "../data/g4470/Lam_4_39", FULL);
    relation* rel_INT1_22_11_10_9_8_7_6_5_4 = new relation(8, false, 22, 259, "rel_INT1_22_11_10_9_8_7_6_5_4", "../data/g4470/INT1_22_38", FULL);
    relation* rel_App_4_1_2_3_4 = new relation(4, true, 4, 273, "rel_App_4_1_2_3_4", "../data/g4470/App_4_37", FULL);
    relation* rel_Free_2_2 = new relation(1, false, 2, 265, "rel_Free_2_2", "../data/g4470/Free_2_36", FULL);
    relation* rel_INT1_22_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22 = new relation(22, true, 22, 259, "rel_INT1_22_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22", "../data/g4470/INT1_22_35", FULL);
    relation* rel_INT2_30_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30 = new relation(30, true, 30, 263, "rel_INT2_30_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30", "../data/g4470/INT2_30_34", FULL);
    relation* rel_Step_24_17 = new relation(1, false, 24, 267, "rel_Step_24_17", "../data/g4470/Step_24_33", FULL);
    relation* rel_AEval_16_8_7_6_5_4_3_2_1 = new relation(8, false, 16, 266, "rel_AEval_16_8_7_6_5_4_3_2_1", "../data/g4470/AEval_16_32", FULL);
    relation* rel_Prog_1_ = new relation(0, false, 1, 260, "rel_Prog_1_", "../data/g4470/Prog_1_31", FULL);
    relation* rel_Free1_3_1_2_3 = new relation(3, true, 3, 261, "rel_Free1_3_1_2_3", "../data/g4470/Free1_3_30", FULL);
    relation* rel_INT0_19_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19 = new relation(19, true, 19, 274, "rel_INT0_19_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19", "../data/g4470/INT0_19_29", FULL);
    relation* rel_Free1_3_2_1 = new relation(2, false, 3, 261, "rel_Free1_3_2_1", "../data/g4470/Free1_3_28", FULL);
    relation* rel_ReachesCfg_8_1 = new relation(1, false, 8, 270, "rel_ReachesCfg_8_1", "../data/g4470/ReachesCfg_8_27", FULL);
    relation* rel_Store_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16 = new relation(16, true, 16, 272, "rel_Store_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16", "../data/g4470/Store_16_26", FULL);
    relation* rel_AE0_9_9_8_7_6_5_4_3_2 = new relation(8, false, 9, 275, "rel_AE0_9_9_8_7_6_5_4_3_2", "../data/g4470/AE0_9_25", FULL);
    relation* rel_App_4_2 = new relation(1, false, 4, 273, "rel_App_4_2", "../data/g4470/App_4_24", FULL);
    relation* rel_INT2_30_11_10_9_8_7_6_5_3 = new relation(8, false, 30, 263, "rel_INT2_30_11_10_9_8_7_6_5_3", "../data/g4470/INT2_30_23", FULL);
    relation* rel_App_4_3 = new relation(1, false, 4, 273, "rel_App_4_3", "../data/g4470/App_4_22", FULL);
    relation* rel_Step_24_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24 = new relation(24, true, 24, 267, "rel_Step_24_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24", "../data/g4470/Step_24_21", FULL);
    relation* rel_AE0_9_1_2_3_4_5_6_7_8_9 = new relation(9, true, 9, 275, "rel_AE0_9_1_2_3_4_5_6_7_8_9", "../data/g4470/AE0_9_20", FULL);
    relation* rel_Lam_4_4 = new relation(1, false, 4, 264, "rel_Lam_4_4", "../data/g4470/Lam_4_19", FULL);
    relation* rel_Lam_4_1 = new relation(1, false, 4, 264, "rel_Lam_4_1", "../data/g4470/Lam_4_18", FULL);
    relation* rel_INT00_11_1_2_3_4_5_6_7_8_9_10_11 = new relation(11, true, 11, 256, "rel_INT00_11_1_2_3_4_5_6_7_8_9_10_11", "../data/g4470/INT00_11_17", FULL);
    relation* rel_FrProp_15_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15 = new relation(15, true, 15, 257, "rel_FrProp_15_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15", "../data/g4470/FrProp_15_16", FULL);
    relation* rel_App_4_4 = new relation(1, false, 4, 273, "rel_App_4_4", "../data/g4470/App_4_15", FULL);
    relation* rel_INT00_11_11_10_9_8_7_6_5_2 = new relation(8, false, 11, 256, "rel_INT00_11_11_10_9_8_7_6_5_2", "../data/g4470/INT00_11_14", FULL);
    relation* rel_Time_7_1_2_3_4_5_6_7 = new relation(7, true, 7, 269, "rel_Time_7_1_2_3_4_5_6_7", "../data/g4470/Time_7_13", FULL);
    relation* rel_Prog_1_1 = new relation(1, true, 1, 260, "rel_Prog_1_1", "../data/g4470/Prog_1_12", FULL);
    relation* rel_INT0_19_12 = new relation(1, false, 19, 274, "rel_INT0_19_12", "../data/g4470/INT0_19_11", FULL);
    relation* rel_AEval_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16 = new relation(16, true, 16, 266, "rel_AEval_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16", "../data/g4470/AEval_16_10", FULL);
    relation* rel_ReachesClo_8_1_2_3_4_5_6_7_8 = new relation(8, true, 8, 276, "rel_ReachesClo_8_1_2_3_4_5_6_7_8", "../data/g4470/ReachesClo_8_9", FULL);
    relation* rel_ReachesCfg_8_1_2_3_4_5_6_7_8 = new relation(8, true, 8, 270, "rel_ReachesCfg_8_1_2_3_4_5_6_7_8", "../data/g4470/ReachesCfg_8_8", FULL);
    relation* rel_Store_16_8_7_6_5_4_3_2_1 = new relation(8, false, 16, 272, "rel_Store_16_8_7_6_5_4_3_2_1", "../data/g4470/Store_16_7", FULL);
    relation* rel_Free0_4_1_2_3_4 = new relation(4, true, 4, 258, "rel_Free0_4_1_2_3_4", "../data/g4470/Free0_4_6", FULL);
    relation* rel_Var_2_1 = new relation(1, false, 2, 268, "rel_Var_2_1", "../data/g4470/Var_2_5", FULL);
    relation* rel_Free_2_1_2 = new relation(2, true, 2, 265, "rel_Free_2_1_2", "../data/g4470/Free_2_4", FULL);
    relation* rel_Var_2_1_2 = new relation(2, true, 2, 268, "rel_Var_2_1_2", "../data/g4470/Var_2_3", FULL);
    relation* rel_APP_35_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35 = new relation(35, true, 35, 262, "rel_APP_35_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35", "../data/g4470/APP_35_2", FULL);
    relation* rel_App_4_1 = new relation(1, false, 4, 273, "rel_App_4_1", "../data/g4470/App_4_1", FULL);
    relation* rel_FrProp_15_14_13_12_11_10_9_8_15 = new relation(8, false, 15, 257, "rel_FrProp_15_14_13_12_11_10_9_8_15", "../data/g4470/FrProp_15_0", FULL);

    RAM* scc4471 = new RAM(false, 1);
    scc4471->add_relation(rel_Lam_4_4, true);
    scc4471->add_relation(rel_Lam_4_1_2_3_4, true);
    scc4471->add_rule(new parallel_acopy(rel_Lam_4_4, rel_Lam_4_1_2_3_4, DELTA, {3, 4, 0, 1, 2}));

    RAM* scc4472 = new RAM(false, 5);
    scc4472->add_relation(rel_App_4_4, true);
    scc4472->add_relation(rel_App_4_1_2_3_4, true);
    scc4472->add_rule(new parallel_acopy(rel_App_4_4, rel_App_4_1_2_3_4, DELTA, {3, 4, 0, 1, 2}));

    RAM* scc4473 = new RAM(false, 9);
    scc4473->add_relation(rel_APP_35_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35, false);
    scc4473->add_relation(rel_ReachesClo_8_1_2_3_4_5_6_7_8, true);
    scc4473->add_rule(new parallel_copy(rel_ReachesClo_8_1_2_3_4_5_6_7_8, rel_APP_35_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35, FULL, {8, 9, 10, 11, 12, 13, 14, 15}));

    RAM* scc4474 = new RAM(false, 13);
    scc4474->add_relation(rel_Prog_1_1, false);
    scc4474->add_relation(rel_Time_7_1_2_3_4_5_6_7, true);
    scc4474->add_rule(new parallel_copy(rel_Time_7_1_2_3_4_5_6_7, rel_Prog_1_1, FULL, {0, 0, 0, 0, 0, 0, 0}));

    RAM* scc4475 = new RAM(false, 3);
    scc4475->add_relation(rel_Var_2_1_2, true);
    scc4475->add_relation(rel_Var_2_1, true);
    scc4475->add_rule(new parallel_acopy(rel_Var_2_1, rel_Var_2_1_2, DELTA, {0, 2, 1}));

    RAM* scc4476 = new RAM(false, 7);
    scc4476->add_relation(rel_ReachesCfg_8_1_2_3_4_5_6_7_8, true);
    scc4476->add_relation(rel_Prog_1_1, false);
    scc4476->add_rule(new parallel_copy(rel_ReachesCfg_8_1_2_3_4_5_6_7_8, rel_Prog_1_1, FULL, {0, 0, 0, 0, 0, 0, 0, 0}));

    RAM* scc4477 = new RAM(false, 11);
    scc4477->add_relation(rel_App_4_1, true);
    scc4477->add_relation(rel_App_4_1_2_3_4, true);
    scc4477->add_rule(new parallel_acopy(rel_App_4_1, rel_App_4_1_2_3_4, DELTA, {0, 4, 1, 2, 3}));

    RAM* scc4478 = new RAM(false, 15);
    scc4478->add_relation(rel_App_4_2, true);
    scc4478->add_relation(rel_App_4_1_2_3_4, true);
    scc4478->add_rule(new parallel_acopy(rel_App_4_2, rel_App_4_1_2_3_4, DELTA, {1, 4, 0, 2, 3}));

    RAM* scc4479 = new RAM(false, 2);
    scc4479->add_relation(rel_APP_35_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35, false);
    scc4479->add_relation(rel_Time_7_1_2_3_4_5_6_7, true);
    scc4479->add_rule(new parallel_copy(rel_Time_7_1_2_3_4_5_6_7, rel_APP_35_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35, FULL, {0, 1, 2, 3, 4, 5, 6}));

    RAM* scc4480 = new RAM(false, 6);
    scc4480->add_relation(rel_Prog_1_1, true);
    scc4480->add_relation(rel_Prog_1_, true);
    scc4480->add_rule(new parallel_acopy(rel_Prog_1_, rel_Prog_1_1, DELTA, {1, 0}));

    RAM* scc4481 = new RAM(true, 10);
    scc4481->add_relation(rel_FrProp_15_14_13_12_11_10_9_8_15, true);
    scc4481->add_relation(rel_App_4_1, false);
    scc4481->add_relation(rel_APP_35_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35, true);
    scc4481->add_relation(rel_Var_2_1, false);
    scc4481->add_relation(rel_Store_16_8_7_6_5_4_3_2_1, true);
    scc4481->add_relation(rel_ReachesCfg_8_1_2_3_4_5_6_7_8, true);
    scc4481->add_relation(rel_AEval_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16, true);
    scc4481->add_relation(rel_INT0_19_12, true);
    scc4481->add_relation(rel_INT00_11_11_10_9_8_7_6_5_2, true);
    scc4481->add_relation(rel_FrProp_15_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15, true);
    scc4481->add_relation(rel_INT00_11_1_2_3_4_5_6_7_8_9_10_11, true);
    scc4481->add_relation(rel_Lam_4_1, false);
    scc4481->add_relation(rel_AE0_9_1_2_3_4_5_6_7_8_9, true);
    scc4481->add_relation(rel_Step_24_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24, true);
    scc4481->add_relation(rel_INT2_30_11_10_9_8_7_6_5_3, true);
    scc4481->add_relation(rel_AE0_9_9_8_7_6_5_4_3_2, true);
    scc4481->add_relation(rel_Store_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16, true);
    scc4481->add_relation(rel_ReachesCfg_8_1, true);
    scc4481->add_relation(rel_INT0_19_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19, true);
    scc4481->add_relation(rel_AEval_16_8_7_6_5_4_3_2_1, true);
    scc4481->add_relation(rel_Step_24_17, true);
    scc4481->add_relation(rel_INT2_30_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30, true);
    scc4481->add_relation(rel_INT1_22_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22, true);
    scc4481->add_relation(rel_Free_2_2, false);
    scc4481->add_relation(rel_INT1_22_11_10_9_8_7_6_5_4, true);
    scc4481->add_rule(new parallel_join(rel_INT00_11_1_2_3_4_5_6_7_8_9_10_11, rel_App_4_1, FULL, rel_ReachesCfg_8_1, DELTA, {0, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12}));
    scc4481->add_rule(new parallel_acopy(rel_AE0_9_9_8_7_6_5_4_3_2, rel_AE0_9_1_2_3_4_5_6_7_8_9, DELTA, {8, 7, 6, 5, 4, 3, 2, 1, 9, 0}));
    scc4481->add_rule(new parallel_acopy(rel_INT0_19_12, rel_INT0_19_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19, DELTA, {11, 19, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18}));
    scc4481->add_rule(new parallel_join(rel_ReachesCfg_8_1_2_3_4_5_6_7_8, rel_App_4_1, FULL, rel_ReachesCfg_8_1, DELTA, {2, 6, 7, 8, 9, 10, 11, 12}));
    scc4481->add_rule(new parallel_copy(rel_Store_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16, rel_APP_35_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35, DELTA, {18, 0, 1, 2, 3, 4, 5, 6, 27, 28, 29, 30, 31, 32, 33, 34}));
    scc4481->add_rule(new parallel_acopy(rel_INT2_30_11_10_9_8_7_6_5_3, rel_INT2_30_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30, DELTA, {10, 9, 8, 7, 6, 5, 4, 2, 30, 0, 1, 3, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29}));
    scc4481->add_rule(new parallel_join(rel_AEval_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16, rel_Lam_4_1, FULL, rel_ReachesCfg_8_1, DELTA, {0, 6, 7, 8, 9, 10, 11, 12, 0, 6, 7, 8, 9, 10, 11, 12}));
    scc4481->add_rule(new parallel_copy(rel_Store_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16, rel_APP_35_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35, DELTA, {17, 0, 1, 2, 3, 4, 5, 6, 19, 20, 21, 22, 23, 24, 25, 26}));
    scc4481->add_rule(new parallel_acopy(rel_AEval_16_8_7_6_5_4_3_2_1, rel_AEval_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16, DELTA, {7, 6, 5, 4, 3, 2, 1, 0, 16, 8, 9, 10, 11, 12, 13, 14, 15}));
    scc4481->add_rule(new parallel_acopy(rel_Step_24_17, rel_Step_24_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24, DELTA, {16, 24, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23}));
    scc4481->add_rule(new parallel_join(rel_AE0_9_1_2_3_4_5_6_7_8_9, rel_Var_2_1, FULL, rel_ReachesCfg_8_1, DELTA, {0, 2, 4, 5, 6, 7, 8, 9, 10}));
    scc4481->add_rule(new parallel_join(rel_Store_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16, rel_FrProp_15_14_13_12_11_10_9_8_15, FULL, rel_Store_16_8_7_6_5_4_3_2_1, DELTA, {7, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24}));
    scc4481->add_rule(new parallel_join(rel_ReachesCfg_8_1_2_3_4_5_6_7_8, rel_App_4_1, FULL, rel_ReachesCfg_8_1, DELTA, {3, 6, 7, 8, 9, 10, 11, 12}));
    scc4481->add_rule(new parallel_acopy(rel_FrProp_15_14_13_12_11_10_9_8_15, rel_FrProp_15_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15, DELTA, {13, 12, 11, 10, 9, 8, 7, 14, 15, 0, 1, 2, 3, 4, 5, 6}));
    scc4481->add_rule(new parallel_join(rel_INT1_22_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22, rel_Lam_4_1, FULL, rel_INT0_19_12, DELTA, {6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 0, 17, 18, 19, 20, 21, 22, 23, 2, 3, 4}));
    scc4481->add_rule(new parallel_acopy(rel_INT1_22_11_10_9_8_7_6_5_4, rel_INT1_22_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22, DELTA, {10, 9, 8, 7, 6, 5, 4, 3, 22, 0, 1, 2, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21}));
    scc4481->add_rule(new parallel_copy(rel_Step_24_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24, rel_APP_35_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35, DELTA, {0, 1, 2, 3, 4, 5, 6, 7, 16, 0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15}));
    scc4481->add_rule(new parallel_join(rel_INT0_19_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19, rel_INT00_11_11_10_9_8_7_6_5_2, DELTA, rel_AEval_16_8_7_6_5_4_3_2_1, FULL, {9, 7, 10, 11, 6, 5, 4, 3, 2, 1, 0, 13, 14, 15, 16, 17, 18, 19, 20}));
    scc4481->add_rule(new parallel_join(rel_INT0_19_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19, rel_INT00_11_11_10_9_8_7_6_5_2, DELTA, rel_AEval_16_8_7_6_5_4_3_2_1, DELTA, {9, 7, 10, 11, 6, 5, 4, 3, 2, 1, 0, 13, 14, 15, 16, 17, 18, 19, 20}));
    scc4481->add_rule(new parallel_join(rel_Store_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16, rel_FrProp_15_14_13_12_11_10_9_8_15, DELTA, rel_Store_16_8_7_6_5_4_3_2_1, FULL, {7, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24}));
    scc4481->add_rule(new parallel_join(rel_Store_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16, rel_FrProp_15_14_13_12_11_10_9_8_15, DELTA, rel_Store_16_8_7_6_5_4_3_2_1, DELTA, {7, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24}));
    scc4481->add_rule(new parallel_acopy(rel_Store_16_8_7_6_5_4_3_2_1, rel_Store_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16, DELTA, {7, 6, 5, 4, 3, 2, 1, 0, 16, 8, 9, 10, 11, 12, 13, 14, 15}));
    scc4481->add_rule(new parallel_acopy(rel_INT00_11_11_10_9_8_7_6_5_2, rel_INT00_11_1_2_3_4_5_6_7_8_9_10_11, DELTA, {10, 9, 8, 7, 6, 5, 4, 1, 11, 0, 2, 3}));
    scc4481->add_rule(new parallel_join(rel_ReachesCfg_8_1_2_3_4_5_6_7_8, rel_App_4_1, FULL, rel_ReachesCfg_8_1, DELTA, {4, 6, 7, 8, 9, 10, 11, 12}));
    scc4481->add_rule(new parallel_join(rel_AEval_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16, rel_AE0_9_9_8_7_6_5_4_3_2, DELTA, rel_Store_16_8_7_6_5_4_3_2_1, FULL, {9, 6, 5, 4, 3, 2, 1, 0, 11, 12, 13, 14, 15, 16, 17, 18}));
    scc4481->add_rule(new parallel_join(rel_AEval_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16, rel_AE0_9_9_8_7_6_5_4_3_2, DELTA, rel_Store_16_8_7_6_5_4_3_2_1, DELTA, {9, 6, 5, 4, 3, 2, 1, 0, 11, 12, 13, 14, 15, 16, 17, 18}));
    scc4481->add_rule(new parallel_join(rel_AEval_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16, rel_AE0_9_9_8_7_6_5_4_3_2, FULL, rel_Store_16_8_7_6_5_4_3_2_1, DELTA, {9, 6, 5, 4, 3, 2, 1, 0, 11, 12, 13, 14, 15, 16, 17, 18}));
    scc4481->add_rule(new parallel_join(rel_INT0_19_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19, rel_INT00_11_11_10_9_8_7_6_5_2, FULL, rel_AEval_16_8_7_6_5_4_3_2_1, DELTA, {9, 7, 10, 11, 6, 5, 4, 3, 2, 1, 0, 13, 14, 15, 16, 17, 18, 19, 20}));
    scc4481->add_rule(new parallel_join(rel_INT2_30_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30, rel_AEval_16_8_7_6_5_4_3_2_1, DELTA, rel_INT1_22_11_10_9_8_7_6_5_4, FULL, {18, 19, 20, 7, 6, 5, 4, 3, 2, 1, 0, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 9, 10, 11, 12, 13, 14, 15, 16}));
    scc4481->add_rule(new parallel_join(rel_INT2_30_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30, rel_AEval_16_8_7_6_5_4_3_2_1, FULL, rel_INT1_22_11_10_9_8_7_6_5_4, DELTA, {18, 19, 20, 7, 6, 5, 4, 3, 2, 1, 0, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 9, 10, 11, 12, 13, 14, 15, 16}));
    scc4481->add_rule(new parallel_join(rel_INT2_30_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30, rel_AEval_16_8_7_6_5_4_3_2_1, DELTA, rel_INT1_22_11_10_9_8_7_6_5_4, DELTA, {18, 19, 20, 7, 6, 5, 4, 3, 2, 1, 0, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 9, 10, 11, 12, 13, 14, 15, 16}));
    scc4481->add_rule(new parallel_acopy(rel_ReachesCfg_8_1, rel_ReachesCfg_8_1_2_3_4_5_6_7_8, DELTA, {0, 8, 1, 2, 3, 4, 5, 6, 7}));
    scc4481->add_rule(new parallel_join(rel_APP_35_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35, rel_AEval_16_8_7_6_5_4_3_2_1, DELTA, rel_INT2_30_11_10_9_8_7_6_5_3, FULL, {18, 6, 5, 4, 3, 2, 1, 0, 21, 22, 23, 24, 25, 26, 27, 28, 31, 29, 30, 9, 10, 11, 12, 13, 14, 15, 16, 32, 33, 34, 35, 36, 37, 38, 39}));
    scc4481->add_rule(new parallel_join(rel_APP_35_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35, rel_AEval_16_8_7_6_5_4_3_2_1, FULL, rel_INT2_30_11_10_9_8_7_6_5_3, DELTA, {18, 6, 5, 4, 3, 2, 1, 0, 21, 22, 23, 24, 25, 26, 27, 28, 31, 29, 30, 9, 10, 11, 12, 13, 14, 15, 16, 32, 33, 34, 35, 36, 37, 38, 39}));
    scc4481->add_rule(new parallel_join(rel_APP_35_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35, rel_AEval_16_8_7_6_5_4_3_2_1, DELTA, rel_INT2_30_11_10_9_8_7_6_5_3, DELTA, {18, 6, 5, 4, 3, 2, 1, 0, 21, 22, 23, 24, 25, 26, 27, 28, 31, 29, 30, 9, 10, 11, 12, 13, 14, 15, 16, 32, 33, 34, 35, 36, 37, 38, 39}));
    scc4481->add_rule(new parallel_join(rel_FrProp_15_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15, rel_Free_2_2, FULL, rel_Step_24_17, DELTA, {13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 2}));
    scc4481->add_rule(new parallel_copy(rel_ReachesCfg_8_1_2_3_4_5_6_7_8, rel_APP_35_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35, DELTA, {16, 0, 1, 2, 3, 4, 5, 6}));

    RAM* scc4482 = new RAM(false, 14);
    scc4482->add_relation(rel_App_4_3, true);
    scc4482->add_relation(rel_App_4_1_2_3_4, true);
    scc4482->add_rule(new parallel_acopy(rel_App_4_3, rel_App_4_1_2_3_4, DELTA, {2, 4, 0, 1, 3}));

    RAM* scc4483 = new RAM(true, 4);
    scc4483->add_relation(rel_Free_2_1_2, true);
    scc4483->add_relation(rel_Free0_4_1_2_3_4, true);
    scc4483->add_relation(rel_App_4_4, false);
    scc4483->add_relation(rel_Lam_4_4, false);
    scc4483->add_relation(rel_App_4_3, false);
    scc4483->add_relation(rel_App_4_2, false);
    scc4483->add_relation(rel_Free1_3_2_1, true);
    scc4483->add_relation(rel_Free1_3_1_2_3, true);
    scc4483->add_relation(rel_Free_2_2, true);
    scc4483->add_relation(rel_Free0_4_2_1, true);
    scc4483->add_rule(new parallel_join(rel_Free_2_1_2, rel_Free_2_2, DELTA, rel_App_4_2, FULL, {2, 4}));
    scc4483->add_rule(new parallel_acopy(rel_Free0_4_2_1, rel_Free0_4_1_2_3_4, DELTA, {1, 0, 4, 2, 3}));
    scc4483->add_rule(new parallel_join(rel_Free_2_1_2, rel_Free_2_2, DELTA, rel_App_4_4, FULL, {2, 4}));
    scc4483->add_rule(new parallel_copy_filter(rel_Free_2_1_2, rel_Free1_3_2_1, DELTA, {1, 3}, [](const u64* const data){ return !(data[0] == data[1]); }));
    scc4483->add_rule(new parallel_copy_filter(rel_Free1_3_1_2_3, rel_Free0_4_2_1, DELTA, {1, 3, 4}, [](const u64* const data){ return !(data[0] == data[1]); }));
    scc4483->add_rule(new parallel_acopy(rel_Free1_3_2_1, rel_Free1_3_1_2_3, DELTA, {1, 0, 3, 2}));
    scc4483->add_rule(new parallel_acopy(rel_Free_2_2, rel_Free_2_1_2, DELTA, {1, 2, 0}));
    scc4483->add_rule(new parallel_join(rel_Free_2_1_2, rel_Free_2_2, DELTA, rel_App_4_3, FULL, {2, 4}));
    scc4483->add_rule(new parallel_join(rel_Free0_4_1_2_3_4, rel_Free_2_2, DELTA, rel_Lam_4_4, FULL, {2, 5, 6, 4}));

    RAM* scc4484 = new RAM(false, 8);
    scc4484->add_relation(rel_Var_2_1_2, false);
    scc4484->add_relation(rel_Free_2_1_2, true);
    scc4484->add_rule(new parallel_copy(rel_Free_2_1_2, rel_Var_2_1_2, FULL, {1, 0}));

    RAM* scc4485 = new RAM(false, 12);
    scc4485->add_relation(rel_Lam_4_1, true);
    scc4485->add_relation(rel_Lam_4_1_2_3_4, true);
    scc4485->add_rule(new parallel_acopy(rel_Lam_4_1, rel_Lam_4_1_2_3_4, DELTA, {0, 4, 1, 2, 3}));

    LIE* lie = new LIE();
    lie->add_relation(rel_Free0_4_2_1);
    lie->add_relation(rel_Lam_4_1_2_3_4);
    lie->add_relation(rel_INT1_22_11_10_9_8_7_6_5_4);
    lie->add_relation(rel_App_4_1_2_3_4);
    lie->add_relation(rel_Free_2_2);
    lie->add_relation(rel_INT1_22_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22);
    lie->add_relation(rel_INT2_30_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30);
    lie->add_relation(rel_Step_24_17);
    lie->add_relation(rel_AEval_16_8_7_6_5_4_3_2_1);
    lie->add_relation(rel_Prog_1_);
    lie->add_relation(rel_Free1_3_1_2_3);
    lie->add_relation(rel_INT0_19_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19);
    lie->add_relation(rel_Free1_3_2_1);
    lie->add_relation(rel_ReachesCfg_8_1);
    lie->add_relation(rel_Store_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16);
    lie->add_relation(rel_AE0_9_9_8_7_6_5_4_3_2);
    lie->add_relation(rel_App_4_2);
    lie->add_relation(rel_INT2_30_11_10_9_8_7_6_5_3);
    lie->add_relation(rel_App_4_3);
    lie->add_relation(rel_Step_24_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24);
    lie->add_relation(rel_AE0_9_1_2_3_4_5_6_7_8_9);
    lie->add_relation(rel_Lam_4_4);
    lie->add_relation(rel_Lam_4_1);
    lie->add_relation(rel_INT00_11_1_2_3_4_5_6_7_8_9_10_11);
    lie->add_relation(rel_FrProp_15_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15);
    lie->add_relation(rel_App_4_4);
    lie->add_relation(rel_INT00_11_11_10_9_8_7_6_5_2);
    lie->add_relation(rel_Time_7_1_2_3_4_5_6_7);
    lie->add_relation(rel_Prog_1_1);
    lie->add_relation(rel_INT0_19_12);
    lie->add_relation(rel_AEval_16_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16);
    lie->add_relation(rel_ReachesClo_8_1_2_3_4_5_6_7_8);
    lie->add_relation(rel_ReachesCfg_8_1_2_3_4_5_6_7_8);
    lie->add_relation(rel_Store_16_8_7_6_5_4_3_2_1);
    lie->add_relation(rel_Free0_4_1_2_3_4);
    lie->add_relation(rel_Var_2_1);
    lie->add_relation(rel_Free_2_1_2);
    lie->add_relation(rel_Var_2_1_2);
    lie->add_relation(rel_APP_35_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35);
    lie->add_relation(rel_App_4_1);
    lie->add_relation(rel_FrProp_15_14_13_12_11_10_9_8_15);
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
    lie->add_scc_dependance(scc4471, scc4483);
    lie->add_scc_dependance(scc4472, scc4483);
    lie->add_scc_dependance(scc4475, scc4481);
    lie->add_scc_dependance(scc4476, scc4481);
    lie->add_scc_dependance(scc4477, scc4481);
    lie->add_scc_dependance(scc4478, scc4483);
    lie->add_scc_dependance(scc4481, scc4479);
    lie->add_scc_dependance(scc4481, scc4473);
    lie->add_scc_dependance(scc4482, scc4483);
    lie->add_scc_dependance(scc4483, scc4481);
    lie->add_scc_dependance(scc4484, scc4483);
    lie->add_scc_dependance(scc4485, scc4481);



    lie->set_comm(mcomm);
    lie->set_batch_size(1);
    lie->execute();

    delete lie;

    mcomm.destroy();
    return 0;
}
