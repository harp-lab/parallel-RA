#include "../src/parallel_RA_inc.h"


int main(int argc, char **argv)
{
    mpi_comm mcomm;
    mcomm.create(argc, argv);



    relation* rel_Free0_4_2_1 = new relation(2, false, 4, 258, "rel_Free0_4_2_1", "../data/g4470/Free0_4_40", FULL);
    relation* rel_AE0_10_10_9_8_7_6_5_4_3_2 = new relation(9, false, 10, 268, "rel_AE0_10_10_9_8_7_6_5_4_3_2", "../data/g4470/AE0_10_39", FULL);
    relation* rel_Lam_4_1_2_3_4 = new relation(4, true, 4, 266, "rel_Lam_4_1_2_3_4", "../data/g4470/Lam_4_38", FULL);
    relation* rel_App_4_1_2_3_4 = new relation(4, true, 4, 274, "rel_App_4_1_2_3_4", "../data/g4470/App_4_37", FULL);
    relation* rel_Free_2_2 = new relation(1, false, 2, 267, "rel_Free_2_2", "../data/g4470/Free_2_36", FULL);
    relation* rel_Step_27_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27 = new relation(27, true, 27, 265, "rel_Step_27_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27", "../data/g4470/Step_27_35", FULL);
    relation* rel_AE0_10_1_2_3_4_5_6_7_8_9_10 = new relation(10, true, 10, 268, "rel_AE0_10_1_2_3_4_5_6_7_8_9_10", "../data/g4470/AE0_10_34", FULL);
    relation* rel_Step_27_19 = new relation(1, false, 27, 265, "rel_Step_27_19", "../data/g4470/Step_27_33", FULL);
    relation* rel_Prog_1_ = new relation(0, false, 1, 259, "rel_Prog_1_", "../data/g4470/Prog_1_32", FULL);
    relation* rel_Free1_3_1_2_3 = new relation(3, true, 3, 260, "rel_Free1_3_1_2_3", "../data/g4470/Free1_3_31", FULL);
    relation* rel_Free1_3_2_1 = new relation(2, false, 3, 260, "rel_Free1_3_2_1", "../data/g4470/Free1_3_30", FULL);
    relation* rel_INT00_12_1_2_3_4_5_6_7_8_9_10_11_12 = new relation(12, true, 12, 271, "rel_INT00_12_1_2_3_4_5_6_7_8_9_10_11_12", "../data/g4470/INT00_12_29", FULL);
    relation* rel_ReachesClo_9_1_2_3_4_5_6_7_8_9 = new relation(9, true, 9, 269, "rel_ReachesClo_9_1_2_3_4_5_6_7_8_9", "../data/g4470/ReachesClo_9_28", FULL);
    relation* rel_INT2_33_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33 = new relation(33, true, 33, 272, "rel_INT2_33_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33", "../data/g4470/INT2_33_27", FULL);
    relation* rel_App_4_2 = new relation(1, false, 4, 274, "rel_App_4_2", "../data/g4470/App_4_26", FULL);
    relation* rel_INT1_24_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24 = new relation(24, true, 24, 264, "rel_INT1_24_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24", "../data/g4470/INT1_24_25", FULL);
    relation* rel_Store_18_9_8_7_6_5_4_3_2_1 = new relation(9, false, 18, 276, "rel_Store_18_9_8_7_6_5_4_3_2_1", "../data/g4470/Store_18_24", FULL);
    relation* rel_App_4_3 = new relation(1, false, 4, 274, "rel_App_4_3", "../data/g4470/App_4_23", FULL);
    relation* rel_FrProp_17_16_15_14_13_12_11_10_9_17 = new relation(9, false, 17, 263, "rel_FrProp_17_16_15_14_13_12_11_10_9_17", "../data/g4470/FrProp_17_22", FULL);
    relation* rel_ReachesCfg_9_1_2_3_4_5_6_7_8_9 = new relation(9, true, 9, 262, "rel_ReachesCfg_9_1_2_3_4_5_6_7_8_9", "../data/g4470/ReachesCfg_9_21", FULL);
    relation* rel_INT1_24_12_11_10_9_8_7_6_5_4 = new relation(9, false, 24, 264, "rel_INT1_24_12_11_10_9_8_7_6_5_4", "../data/g4470/INT1_24_20", FULL);
    relation* rel_Time_8_1_2_3_4_5_6_7_8 = new relation(8, true, 8, 257, "rel_Time_8_1_2_3_4_5_6_7_8", "../data/g4470/Time_8_19", FULL);
    relation* rel_Lam_4_4 = new relation(1, false, 4, 266, "rel_Lam_4_4", "../data/g4470/Lam_4_18", FULL);
    relation* rel_Lam_4_1 = new relation(1, false, 4, 266, "rel_Lam_4_1", "../data/g4470/Lam_4_17", FULL);
    relation* rel_AEval_18_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18 = new relation(18, true, 18, 256, "rel_AEval_18_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18", "../data/g4470/AEval_18_16", FULL);
    relation* rel_INT2_33_12_11_10_9_8_7_6_5_3 = new relation(9, false, 33, 272, "rel_INT2_33_12_11_10_9_8_7_6_5_3", "../data/g4470/INT2_33_15", FULL);
    relation* rel_App_4_4 = new relation(1, false, 4, 274, "rel_App_4_4", "../data/g4470/App_4_14", FULL);
    relation* rel_Prog_1_1 = new relation(1, true, 1, 259, "rel_Prog_1_1", "../data/g4470/Prog_1_13", FULL);
    relation* rel_INT0_21_13 = new relation(1, false, 21, 275, "rel_INT0_21_13", "../data/g4470/INT0_21_12", FULL);
    relation* rel_ReachesCfg_9_1 = new relation(1, false, 9, 262, "rel_ReachesCfg_9_1", "../data/g4470/ReachesCfg_9_11", FULL);
    relation* rel_AEval_18_9_8_7_6_5_4_3_2_1 = new relation(9, false, 18, 256, "rel_AEval_18_9_8_7_6_5_4_3_2_1", "../data/g4470/AEval_18_10", FULL);
    relation* rel_APP_39_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35_36_37_38_39 = new relation(39, true, 39, 261, "rel_APP_39_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35_36_37_38_39", "../data/g4470/APP_39_9", FULL);
    relation* rel_Free0_4_1_2_3_4 = new relation(4, true, 4, 258, "rel_Free0_4_1_2_3_4", "../data/g4470/Free0_4_8", FULL);
    relation* rel_Var_2_1 = new relation(1, false, 2, 270, "rel_Var_2_1", "../data/g4470/Var_2_7", FULL);
    relation* rel_Free_2_1_2 = new relation(2, true, 2, 267, "rel_Free_2_1_2", "../data/g4470/Free_2_6", FULL);
    relation* rel_FrProp_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17 = new relation(17, true, 17, 263, "rel_FrProp_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17", "../data/g4470/FrProp_17_5", FULL);
    relation* rel_Var_2_1_2 = new relation(2, true, 2, 270, "rel_Var_2_1_2", "../data/g4470/Var_2_4", FULL);
    relation* rel_INT00_12_12_11_10_9_8_7_6_5_2 = new relation(9, false, 12, 271, "rel_INT00_12_12_11_10_9_8_7_6_5_2", "../data/g4470/INT00_12_3", FULL);
    relation* rel_App_4_1 = new relation(1, false, 4, 274, "rel_App_4_1", "../data/g4470/App_4_2", FULL);
    relation* rel_INT0_21_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21 = new relation(21, true, 21, 275, "rel_INT0_21_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21", "../data/g4470/INT0_21_1", FULL);
    relation* rel_Store_18_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18 = new relation(18, true, 18, 276, "rel_Store_18_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18", "../data/g4470/Store_18_0", FULL);

    RAM* scc4471 = new RAM(false, 1);
    scc4471->add_relation(rel_APP_39_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35_36_37_38_39, false);
    scc4471->add_relation(rel_Time_8_1_2_3_4_5_6_7_8, true);
    scc4471->add_rule(new parallel_copy(rel_Time_8_1_2_3_4_5_6_7_8, rel_APP_39_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35_36_37_38_39, FULL, {0, 1, 2, 3, 4, 5, 6, 7}));

    RAM* scc4472 = new RAM(false, 5);
    scc4472->add_relation(rel_Var_2_1_2, true);
    scc4472->add_relation(rel_Var_2_1, true);
    scc4472->add_rule(new parallel_acopy(rel_Var_2_1, rel_Var_2_1_2, DELTA, {0, 2, 1}));

    RAM* scc4473 = new RAM(false, 9);
    scc4473->add_relation(rel_Prog_1_1, false);
    scc4473->add_relation(rel_Time_8_1_2_3_4_5_6_7_8, true);
    scc4473->add_rule(new parallel_copy(rel_Time_8_1_2_3_4_5_6_7_8, rel_Prog_1_1, FULL, {0, 0, 0, 0, 0, 0, 0, 0}));

    RAM* scc4474 = new RAM(false, 13);
    scc4474->add_relation(rel_Var_2_1_2, false);
    scc4474->add_relation(rel_Free_2_1_2, true);
    scc4474->add_rule(new parallel_copy(rel_Free_2_1_2, rel_Var_2_1_2, FULL, {1, 0}));

    RAM* scc4475 = new RAM(false, 3);
    scc4475->add_relation(rel_Prog_1_1, false);
    scc4475->add_relation(rel_ReachesCfg_9_1_2_3_4_5_6_7_8_9, true);
    scc4475->add_rule(new parallel_copy(rel_ReachesCfg_9_1_2_3_4_5_6_7_8_9, rel_Prog_1_1, FULL, {0, 0, 0, 0, 0, 0, 0, 0, 0}));

    RAM* scc4476 = new RAM(false, 7);
    scc4476->add_relation(rel_Prog_1_1, true);
    scc4476->add_relation(rel_Prog_1_, true);
    scc4476->add_rule(new parallel_acopy(rel_Prog_1_, rel_Prog_1_1, DELTA, {1, 0}));

    RAM* scc4477 = new RAM(false, 11);
    scc4477->add_relation(rel_Lam_4_1, true);
    scc4477->add_relation(rel_Lam_4_1_2_3_4, true);
    scc4477->add_rule(new parallel_acopy(rel_Lam_4_1, rel_Lam_4_1_2_3_4, DELTA, {0, 4, 1, 2, 3}));

    RAM* scc4478 = new RAM(false, 15);
    scc4478->add_relation(rel_App_4_2, true);
    scc4478->add_relation(rel_App_4_1_2_3_4, true);
    scc4478->add_rule(new parallel_acopy(rel_App_4_2, rel_App_4_1_2_3_4, DELTA, {1, 4, 0, 2, 3}));

    RAM* scc4479 = new RAM(false, 2);
    scc4479->add_relation(rel_Lam_4_4, true);
    scc4479->add_relation(rel_Lam_4_1_2_3_4, true);
    scc4479->add_rule(new parallel_acopy(rel_Lam_4_4, rel_Lam_4_1_2_3_4, DELTA, {3, 4, 0, 1, 2}));

    RAM* scc4480 = new RAM(false, 6);
    scc4480->add_relation(rel_App_4_4, true);
    scc4480->add_relation(rel_App_4_1_2_3_4, true);
    scc4480->add_rule(new parallel_acopy(rel_App_4_4, rel_App_4_1_2_3_4, DELTA, {3, 4, 0, 1, 2}));

    RAM* scc4481 = new RAM(false, 10);
    scc4481->add_relation(rel_App_4_1, true);
    scc4481->add_relation(rel_App_4_1_2_3_4, true);
    scc4481->add_rule(new parallel_acopy(rel_App_4_1, rel_App_4_1_2_3_4, DELTA, {0, 4, 1, 2, 3}));

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
    scc4483->add_rule(new parallel_acopy(rel_Free0_4_2_1, rel_Free0_4_1_2_3_4, DELTA, {1, 0, 4, 2, 3}));
    scc4483->add_rule(new parallel_copy_filter(rel_Free_2_1_2, rel_Free1_3_2_1, DELTA, {1, 3}, [](const u64* const data){ return !(data[0] == data[1]); }));
    scc4483->add_rule(new parallel_join(rel_Free_2_1_2, rel_Free_2_2, DELTA, rel_App_4_3, FULL, {2, 4}));
    scc4483->add_rule(new parallel_join(rel_Free0_4_1_2_3_4, rel_Free_2_2, DELTA, rel_Lam_4_4, FULL, {2, 5, 6, 4}));
    scc4483->add_rule(new parallel_copy_filter(rel_Free1_3_1_2_3, rel_Free0_4_2_1, DELTA, {1, 3, 4}, [](const u64* const data){ return !(data[0] == data[1]); }));
    scc4483->add_rule(new parallel_join(rel_Free_2_1_2, rel_Free_2_2, DELTA, rel_App_4_4, FULL, {2, 4}));
    scc4483->add_rule(new parallel_acopy(rel_Free1_3_2_1, rel_Free1_3_1_2_3, DELTA, {1, 0, 3, 2}));
    scc4483->add_rule(new parallel_join(rel_Free_2_1_2, rel_Free_2_2, DELTA, rel_App_4_2, FULL, {2, 4}));
    scc4483->add_rule(new parallel_acopy(rel_Free_2_2, rel_Free_2_1_2, DELTA, {1, 2, 0}));

    RAM* scc4484 = new RAM(true, 8);
    scc4484->add_relation(rel_Store_18_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18, true);
    scc4484->add_relation(rel_INT0_21_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21, true);
    scc4484->add_relation(rel_App_4_1, false);
    scc4484->add_relation(rel_INT00_12_12_11_10_9_8_7_6_5_2, true);
    scc4484->add_relation(rel_FrProp_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17, true);
    scc4484->add_relation(rel_Var_2_1, false);
    scc4484->add_relation(rel_APP_39_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35_36_37_38_39, true);
    scc4484->add_relation(rel_AEval_18_9_8_7_6_5_4_3_2_1, true);
    scc4484->add_relation(rel_ReachesCfg_9_1, true);
    scc4484->add_relation(rel_INT0_21_13, true);
    scc4484->add_relation(rel_INT2_33_12_11_10_9_8_7_6_5_3, true);
    scc4484->add_relation(rel_AEval_18_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18, true);
    scc4484->add_relation(rel_Lam_4_1, false);
    scc4484->add_relation(rel_INT1_24_12_11_10_9_8_7_6_5_4, true);
    scc4484->add_relation(rel_ReachesCfg_9_1_2_3_4_5_6_7_8_9, true);
    scc4484->add_relation(rel_FrProp_17_16_15_14_13_12_11_10_9_17, true);
    scc4484->add_relation(rel_Store_18_9_8_7_6_5_4_3_2_1, true);
    scc4484->add_relation(rel_INT1_24_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24, true);
    scc4484->add_relation(rel_INT2_33_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33, true);
    scc4484->add_relation(rel_INT00_12_1_2_3_4_5_6_7_8_9_10_11_12, true);
    scc4484->add_relation(rel_Step_27_19, true);
    scc4484->add_relation(rel_AE0_10_1_2_3_4_5_6_7_8_9_10, true);
    scc4484->add_relation(rel_Step_27_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27, true);
    scc4484->add_relation(rel_Free_2_2, false);
    scc4484->add_relation(rel_AE0_10_10_9_8_7_6_5_4_3_2, true);
    scc4484->add_rule(new parallel_acopy(rel_AE0_10_10_9_8_7_6_5_4_3_2, rel_AE0_10_1_2_3_4_5_6_7_8_9_10, DELTA, {9, 8, 7, 6, 5, 4, 3, 2, 1, 10, 0}));
    scc4484->add_rule(new parallel_acopy(rel_INT0_21_13, rel_INT0_21_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21, DELTA, {12, 21, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 20}));
    scc4484->add_rule(new parallel_join(rel_INT00_12_1_2_3_4_5_6_7_8_9_10_11_12, rel_App_4_1, FULL, rel_ReachesCfg_9_1, DELTA, {0, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13}));
    scc4484->add_rule(new parallel_join(rel_ReachesCfg_9_1_2_3_4_5_6_7_8_9, rel_App_4_1, FULL, rel_ReachesCfg_9_1, DELTA, {4, 6, 7, 8, 9, 10, 11, 12, 13}));
    scc4484->add_rule(new parallel_join(rel_INT2_33_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33, rel_AEval_18_9_8_7_6_5_4_3_2_1, DELTA, rel_INT1_24_12_11_10_9_8_7_6_5_4, FULL, {20, 21, 22, 8, 7, 6, 5, 4, 3, 2, 1, 0, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 10, 11, 12, 13, 14, 15, 16, 17, 18}));
    scc4484->add_rule(new parallel_join(rel_INT2_33_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33, rel_AEval_18_9_8_7_6_5_4_3_2_1, FULL, rel_INT1_24_12_11_10_9_8_7_6_5_4, DELTA, {20, 21, 22, 8, 7, 6, 5, 4, 3, 2, 1, 0, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 10, 11, 12, 13, 14, 15, 16, 17, 18}));
    scc4484->add_rule(new parallel_join(rel_INT2_33_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33, rel_AEval_18_9_8_7_6_5_4_3_2_1, DELTA, rel_INT1_24_12_11_10_9_8_7_6_5_4, DELTA, {20, 21, 22, 8, 7, 6, 5, 4, 3, 2, 1, 0, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 10, 11, 12, 13, 14, 15, 16, 17, 18}));
    scc4484->add_rule(new parallel_acopy(rel_AEval_18_9_8_7_6_5_4_3_2_1, rel_AEval_18_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18, DELTA, {8, 7, 6, 5, 4, 3, 2, 1, 0, 18, 9, 10, 11, 12, 13, 14, 15, 16, 17}));
    scc4484->add_rule(new parallel_join(rel_INT0_21_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21, rel_INT00_12_12_11_10_9_8_7_6_5_2, FULL, rel_AEval_18_9_8_7_6_5_4_3_2_1, DELTA, {10, 8, 11, 12, 7, 6, 5, 4, 3, 2, 1, 0, 14, 15, 16, 17, 18, 19, 20, 21, 22}));
    scc4484->add_rule(new parallel_join(rel_AEval_18_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18, rel_AE0_10_10_9_8_7_6_5_4_3_2, DELTA, rel_Store_18_9_8_7_6_5_4_3_2_1, FULL, {10, 7, 6, 5, 4, 3, 2, 1, 0, 12, 13, 14, 15, 16, 17, 18, 19, 20}));
    scc4484->add_rule(new parallel_join(rel_AEval_18_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18, rel_AE0_10_10_9_8_7_6_5_4_3_2, DELTA, rel_Store_18_9_8_7_6_5_4_3_2_1, DELTA, {10, 7, 6, 5, 4, 3, 2, 1, 0, 12, 13, 14, 15, 16, 17, 18, 19, 20}));
    scc4484->add_rule(new parallel_join(rel_ReachesCfg_9_1_2_3_4_5_6_7_8_9, rel_App_4_1, FULL, rel_ReachesCfg_9_1, DELTA, {3, 6, 7, 8, 9, 10, 11, 12, 13}));
    scc4484->add_rule(new parallel_copy(rel_Store_18_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18, rel_APP_39_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35_36_37_38_39, DELTA, {19, 0, 1, 2, 3, 4, 5, 6, 7, 21, 22, 23, 24, 25, 26, 27, 28, 29}));
    scc4484->add_rule(new parallel_join(rel_AEval_18_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18, rel_AE0_10_10_9_8_7_6_5_4_3_2, FULL, rel_Store_18_9_8_7_6_5_4_3_2_1, DELTA, {10, 7, 6, 5, 4, 3, 2, 1, 0, 12, 13, 14, 15, 16, 17, 18, 19, 20}));
    scc4484->add_rule(new parallel_join(rel_Store_18_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18, rel_FrProp_17_16_15_14_13_12_11_10_9_17, DELTA, rel_Store_18_9_8_7_6_5_4_3_2_1, FULL, {8, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27}));
    scc4484->add_rule(new parallel_join(rel_Store_18_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18, rel_FrProp_17_16_15_14_13_12_11_10_9_17, DELTA, rel_Store_18_9_8_7_6_5_4_3_2_1, DELTA, {8, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27}));
    scc4484->add_rule(new parallel_join(rel_INT1_24_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24, rel_Lam_4_1, FULL, rel_INT0_21_13, DELTA, {6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 0, 18, 19, 20, 21, 22, 23, 24, 25, 2, 3, 4}));
    scc4484->add_rule(new parallel_acopy(rel_FrProp_17_16_15_14_13_12_11_10_9_17, rel_FrProp_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17, DELTA, {15, 14, 13, 12, 11, 10, 9, 8, 16, 17, 0, 1, 2, 3, 4, 5, 6, 7}));
    scc4484->add_rule(new parallel_join(rel_APP_39_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35_36_37_38_39, rel_AEval_18_9_8_7_6_5_4_3_2_1, DELTA, rel_INT2_33_12_11_10_9_8_7_6_5_3, FULL, {20, 7, 6, 5, 4, 3, 2, 1, 0, 23, 24, 25, 26, 27, 28, 29, 30, 31, 34, 32, 33, 10, 11, 12, 13, 14, 15, 16, 17, 18, 35, 36, 37, 38, 39, 40, 41, 42, 43}));
    scc4484->add_rule(new parallel_join(rel_APP_39_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35_36_37_38_39, rel_AEval_18_9_8_7_6_5_4_3_2_1, FULL, rel_INT2_33_12_11_10_9_8_7_6_5_3, DELTA, {20, 7, 6, 5, 4, 3, 2, 1, 0, 23, 24, 25, 26, 27, 28, 29, 30, 31, 34, 32, 33, 10, 11, 12, 13, 14, 15, 16, 17, 18, 35, 36, 37, 38, 39, 40, 41, 42, 43}));
    scc4484->add_rule(new parallel_join(rel_APP_39_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35_36_37_38_39, rel_AEval_18_9_8_7_6_5_4_3_2_1, DELTA, rel_INT2_33_12_11_10_9_8_7_6_5_3, DELTA, {20, 7, 6, 5, 4, 3, 2, 1, 0, 23, 24, 25, 26, 27, 28, 29, 30, 31, 34, 32, 33, 10, 11, 12, 13, 14, 15, 16, 17, 18, 35, 36, 37, 38, 39, 40, 41, 42, 43}));
    scc4484->add_rule(new parallel_acopy(rel_INT00_12_12_11_10_9_8_7_6_5_2, rel_INT00_12_1_2_3_4_5_6_7_8_9_10_11_12, DELTA, {11, 10, 9, 8, 7, 6, 5, 4, 1, 12, 0, 2, 3}));
    scc4484->add_rule(new parallel_join(rel_AEval_18_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18, rel_Lam_4_1, FULL, rel_ReachesCfg_9_1, DELTA, {0, 6, 7, 8, 9, 10, 11, 12, 13, 0, 6, 7, 8, 9, 10, 11, 12, 13}));
    scc4484->add_rule(new parallel_join(rel_FrProp_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17, rel_Free_2_2, FULL, rel_Step_27_19, DELTA, {14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 2}));
    scc4484->add_rule(new parallel_join(rel_AE0_10_1_2_3_4_5_6_7_8_9_10, rel_Var_2_1, FULL, rel_ReachesCfg_9_1, DELTA, {0, 2, 4, 5, 6, 7, 8, 9, 10, 11}));
    scc4484->add_rule(new parallel_copy(rel_Store_18_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18, rel_APP_39_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35_36_37_38_39, DELTA, {20, 0, 1, 2, 3, 4, 5, 6, 7, 30, 31, 32, 33, 34, 35, 36, 37, 38}));
    scc4484->add_rule(new parallel_acopy(rel_ReachesCfg_9_1, rel_ReachesCfg_9_1_2_3_4_5_6_7_8_9, DELTA, {0, 9, 1, 2, 3, 4, 5, 6, 7, 8}));
    scc4484->add_rule(new parallel_join(rel_INT0_21_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21, rel_INT00_12_12_11_10_9_8_7_6_5_2, DELTA, rel_AEval_18_9_8_7_6_5_4_3_2_1, FULL, {10, 8, 11, 12, 7, 6, 5, 4, 3, 2, 1, 0, 14, 15, 16, 17, 18, 19, 20, 21, 22}));
    scc4484->add_rule(new parallel_join(rel_INT0_21_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21, rel_INT00_12_12_11_10_9_8_7_6_5_2, DELTA, rel_AEval_18_9_8_7_6_5_4_3_2_1, DELTA, {10, 8, 11, 12, 7, 6, 5, 4, 3, 2, 1, 0, 14, 15, 16, 17, 18, 19, 20, 21, 22}));
    scc4484->add_rule(new parallel_join(rel_Store_18_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18, rel_FrProp_17_16_15_14_13_12_11_10_9_17, FULL, rel_Store_18_9_8_7_6_5_4_3_2_1, DELTA, {8, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27}));
    scc4484->add_rule(new parallel_acopy(rel_INT2_33_12_11_10_9_8_7_6_5_3, rel_INT2_33_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33, DELTA, {11, 10, 9, 8, 7, 6, 5, 4, 2, 33, 0, 1, 3, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32}));
    scc4484->add_rule(new parallel_copy(rel_Step_27_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27, rel_APP_39_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35_36_37_38_39, DELTA, {0, 1, 2, 3, 4, 5, 6, 7, 8, 18, 0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17}));
    scc4484->add_rule(new parallel_acopy(rel_Store_18_9_8_7_6_5_4_3_2_1, rel_Store_18_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18, DELTA, {8, 7, 6, 5, 4, 3, 2, 1, 0, 18, 9, 10, 11, 12, 13, 14, 15, 16, 17}));
    scc4484->add_rule(new parallel_copy(rel_ReachesCfg_9_1_2_3_4_5_6_7_8_9, rel_APP_39_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35_36_37_38_39, DELTA, {18, 0, 1, 2, 3, 4, 5, 6, 7}));
    scc4484->add_rule(new parallel_acopy(rel_Step_27_19, rel_Step_27_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27, DELTA, {18, 27, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26}));
    scc4484->add_rule(new parallel_join(rel_ReachesCfg_9_1_2_3_4_5_6_7_8_9, rel_App_4_1, FULL, rel_ReachesCfg_9_1, DELTA, {2, 6, 7, 8, 9, 10, 11, 12, 13}));
    scc4484->add_rule(new parallel_acopy(rel_INT1_24_12_11_10_9_8_7_6_5_4, rel_INT1_24_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24, DELTA, {11, 10, 9, 8, 7, 6, 5, 4, 3, 24, 0, 1, 2, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23}));

    RAM* scc4485 = new RAM(false, 12);
    scc4485->add_relation(rel_APP_39_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35_36_37_38_39, false);
    scc4485->add_relation(rel_ReachesClo_9_1_2_3_4_5_6_7_8_9, true);
    scc4485->add_rule(new parallel_copy(rel_ReachesClo_9_1_2_3_4_5_6_7_8_9, rel_APP_39_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35_36_37_38_39, FULL, {9, 10, 11, 12, 13, 14, 15, 16, 17}));

    LIE* lie = new LIE();
    lie->add_relation(rel_Free0_4_2_1);
    lie->add_relation(rel_AE0_10_10_9_8_7_6_5_4_3_2);
    lie->add_relation(rel_Lam_4_1_2_3_4);
    lie->add_relation(rel_App_4_1_2_3_4);
    lie->add_relation(rel_Free_2_2);
    lie->add_relation(rel_Step_27_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27);
    lie->add_relation(rel_AE0_10_1_2_3_4_5_6_7_8_9_10);
    lie->add_relation(rel_Step_27_19);
    lie->add_relation(rel_Prog_1_);
    lie->add_relation(rel_Free1_3_1_2_3);
    lie->add_relation(rel_Free1_3_2_1);
    lie->add_relation(rel_INT00_12_1_2_3_4_5_6_7_8_9_10_11_12);
    lie->add_relation(rel_ReachesClo_9_1_2_3_4_5_6_7_8_9);
    lie->add_relation(rel_INT2_33_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33);
    lie->add_relation(rel_App_4_2);
    lie->add_relation(rel_INT1_24_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24);
    lie->add_relation(rel_Store_18_9_8_7_6_5_4_3_2_1);
    lie->add_relation(rel_App_4_3);
    lie->add_relation(rel_FrProp_17_16_15_14_13_12_11_10_9_17);
    lie->add_relation(rel_ReachesCfg_9_1_2_3_4_5_6_7_8_9);
    lie->add_relation(rel_INT1_24_12_11_10_9_8_7_6_5_4);
    lie->add_relation(rel_Time_8_1_2_3_4_5_6_7_8);
    lie->add_relation(rel_Lam_4_4);
    lie->add_relation(rel_Lam_4_1);
    lie->add_relation(rel_AEval_18_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18);
    lie->add_relation(rel_INT2_33_12_11_10_9_8_7_6_5_3);
    lie->add_relation(rel_App_4_4);
    lie->add_relation(rel_Prog_1_1);
    lie->add_relation(rel_INT0_21_13);
    lie->add_relation(rel_ReachesCfg_9_1);
    lie->add_relation(rel_AEval_18_9_8_7_6_5_4_3_2_1);
    lie->add_relation(rel_APP_39_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35_36_37_38_39);
    lie->add_relation(rel_Free0_4_1_2_3_4);
    lie->add_relation(rel_Var_2_1);
    lie->add_relation(rel_Free_2_1_2);
    lie->add_relation(rel_FrProp_17_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17);
    lie->add_relation(rel_Var_2_1_2);
    lie->add_relation(rel_INT00_12_12_11_10_9_8_7_6_5_2);
    lie->add_relation(rel_App_4_1);
    lie->add_relation(rel_INT0_21_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21);
    lie->add_relation(rel_Store_18_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18);
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
    lie->add_scc_dependance(scc4472, scc4484);
    lie->add_scc_dependance(scc4474, scc4483);
    lie->add_scc_dependance(scc4475, scc4484);
    lie->add_scc_dependance(scc4477, scc4484);
    lie->add_scc_dependance(scc4478, scc4483);
    lie->add_scc_dependance(scc4479, scc4483);
    lie->add_scc_dependance(scc4480, scc4483);
    lie->add_scc_dependance(scc4481, scc4484);
    lie->add_scc_dependance(scc4482, scc4483);
    lie->add_scc_dependance(scc4483, scc4484);
    lie->add_scc_dependance(scc4484, scc4485);
    lie->add_scc_dependance(scc4484, scc4471);








    lie->set_debug_output_filename(argv[1]);
    lie->set_comm(mcomm);
    lie->set_batch_size(100);
    lie->execute();

    delete lie;

    mcomm.destroy();
    return 0;
}
