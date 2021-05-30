#include "../src/parallel_RA_inc.h"




int main(int argc, char **argv)
{
    mpi_comm mcomm;
    mcomm.create(argc, argv);


    {
    relation* T1 = new relation(1, true, 2, 257, "T", "../data/g5955/path_2_1_2", FULL);
    relation* G1 = new relation(1, true, 2, 256, "G", argv[1], FULL);
    RAM* join1 = new RAM(true, 1);
    join1->add_relation(G1, false);
    join1->add_relation(T1, true);
    join1->add_rule(new parallel_join(T1, T1, DELTA, G1, FULL, {3, 1}));
    RAM* copy1 = new RAM(false, 0);
    copy1->add_relation(G1, false);
    copy1->add_relation(T1, true);
    copy1->add_rule(new parallel_copy(T1, G1, FULL, {1, 0}));
    LIE* lie1 = new LIE();
    lie1->add_relation(T1);
    lie1->add_relation(G1);
    lie1->add_scc(join1);
    lie1->add_scc(copy1);
    lie1->add_scc_dependance(copy1, join1);
    //lie1->set_sloav_mode(atoi(argv[3]));
    lie1->set_sloav_mode(0);
    lie1->enable_all_to_all_dump();
    lie1->set_output_dir(argv[2]);
    lie1->enable_IO();
    lie1->set_name(argv[0]);
    lie1->set_comm(mcomm);
    lie1->set_batch_size(10);
    lie1->execute();
    lie1->print_all_relation_size();
    delete lie1;
    }


    {
    relation* T2 = new relation(1, true, 2, 257, "T", "../data/g5955/path_2_1_2", FULL);
    relation* G2 = new relation(1, true, 2, 256, "G", argv[1], FULL);
    RAM* join2 = new RAM(true, 1);
    join2->add_relation(G2, false);
    join2->add_relation(T2, true);
    join2->add_rule(new parallel_join(T2, T2, DELTA, G2, FULL, {3, 1}));
    RAM* copy2 = new RAM(false, 0);
    copy2->add_relation(G2, false);
    copy2->add_relation(T2, true);
    copy2->add_rule(new parallel_copy(T2, G2, FULL, {1, 0}));
    LIE* lie2 = new LIE();
    lie2->add_relation(T2);
    lie2->add_relation(G2);
    lie2->add_scc(join2);
    lie2->add_scc(copy2);
    lie2->add_scc_dependance(copy2, join2);
    //lie2->set_sloav_mode(atoi(argv[5]));
    lie2->set_sloav_mode(1);
    lie2->enable_all_to_all_dump();
    lie2->set_output_dir(argv[3]);
    //lie2->set_output_dir(argv[4]);
    lie2->enable_IO();
    lie2->set_name(argv[0]);
    lie2->set_comm(mcomm);
    lie2->set_batch_size(10);
    lie2->execute();
    lie2->print_all_relation_size();
    delete lie2;
    }


    {
    relation* T1 = new relation(1, true, 2, 257, "T", "../data/g5955/path_2_1_2", FULL);
    relation* G1 = new relation(1, true, 2, 256, "G", argv[1], FULL);
    RAM* join1 = new RAM(true, 1);
    join1->add_relation(G1, false);
    join1->add_relation(T1, true);
    join1->add_rule(new parallel_join(T1, T1, DELTA, G1, FULL, {3, 1}));
    RAM* copy1 = new RAM(false, 0);
    copy1->add_relation(G1, false);
    copy1->add_relation(T1, true);
    copy1->add_rule(new parallel_copy(T1, G1, FULL, {1, 0}));
    LIE* lie1 = new LIE();
    lie1->add_relation(T1);
    lie1->add_relation(G1);
    lie1->add_scc(join1);
    lie1->add_scc(copy1);
    lie1->add_scc_dependance(copy1, join1);
    //lie1->set_sloav_mode(atoi(argv[3]));
    lie1->set_sloav_mode(0);
    lie1->enable_all_to_all_dump();
    lie1->set_output_dir(argv[2]);
    lie1->enable_IO();
    lie1->set_name(argv[0]);
    lie1->set_comm(mcomm);
    lie1->set_batch_size(10);
    lie1->execute();
    lie1->print_all_relation_size();
    delete lie1;
    }


    {
    relation* T2 = new relation(1, true, 2, 257, "T", "../data/g5955/path_2_1_2", FULL);
    relation* G2 = new relation(1, true, 2, 256, "G", argv[1], FULL);
    RAM* join2 = new RAM(true, 1);
    join2->add_relation(G2, false);
    join2->add_relation(T2, true);
    join2->add_rule(new parallel_join(T2, T2, DELTA, G2, FULL, {3, 1}));
    RAM* copy2 = new RAM(false, 0);
    copy2->add_relation(G2, false);
    copy2->add_relation(T2, true);
    copy2->add_rule(new parallel_copy(T2, G2, FULL, {1, 0}));
    LIE* lie2 = new LIE();
    lie2->add_relation(T2);
    lie2->add_relation(G2);
    lie2->add_scc(join2);
    lie2->add_scc(copy2);
    lie2->add_scc_dependance(copy2, join2);
    //lie2->set_sloav_mode(atoi(argv[5]));
    lie2->set_sloav_mode(1);
    lie2->enable_all_to_all_dump();
    lie2->set_output_dir(argv[3]);
    //lie2->set_output_dir(argv[4]);
    lie2->enable_IO();
    lie2->set_name(argv[0]);
    lie2->set_comm(mcomm);
    lie2->set_batch_size(10);
    lie2->execute();
    lie2->print_all_relation_size();
    delete lie2;
    }

    {
    relation* T1 = new relation(1, true, 2, 257, "T", "../data/g5955/path_2_1_2", FULL);
    relation* G1 = new relation(1, true, 2, 256, "G", argv[1], FULL);
    RAM* join1 = new RAM(true, 1);
    join1->add_relation(G1, false);
    join1->add_relation(T1, true);
    join1->add_rule(new parallel_join(T1, T1, DELTA, G1, FULL, {3, 1}));
    RAM* copy1 = new RAM(false, 0);
    copy1->add_relation(G1, false);
    copy1->add_relation(T1, true);
    copy1->add_rule(new parallel_copy(T1, G1, FULL, {1, 0}));
    LIE* lie1 = new LIE();
    lie1->add_relation(T1);
    lie1->add_relation(G1);
    lie1->add_scc(join1);
    lie1->add_scc(copy1);
    lie1->add_scc_dependance(copy1, join1);
    //lie1->set_sloav_mode(atoi(argv[3]));
    lie1->set_sloav_mode(0);
    lie1->enable_all_to_all_dump();
    lie1->set_output_dir(argv[2]);
    lie1->enable_IO();
    lie1->set_name(argv[0]);
    lie1->set_comm(mcomm);
    lie1->set_batch_size(10);
    lie1->execute();
    lie1->print_all_relation_size();
    delete lie1;
    }


    {
    relation* T2 = new relation(1, true, 2, 257, "T", "../data/g5955/path_2_1_2", FULL);
    relation* G2 = new relation(1, true, 2, 256, "G", argv[1], FULL);
    RAM* join2 = new RAM(true, 1);
    join2->add_relation(G2, false);
    join2->add_relation(T2, true);
    join2->add_rule(new parallel_join(T2, T2, DELTA, G2, FULL, {3, 1}));
    RAM* copy2 = new RAM(false, 0);
    copy2->add_relation(G2, false);
    copy2->add_relation(T2, true);
    copy2->add_rule(new parallel_copy(T2, G2, FULL, {1, 0}));
    LIE* lie2 = new LIE();
    lie2->add_relation(T2);
    lie2->add_relation(G2);
    lie2->add_scc(join2);
    lie2->add_scc(copy2);
    lie2->add_scc_dependance(copy2, join2);
    //lie2->set_sloav_mode(atoi(argv[5]));
    lie2->set_sloav_mode(1);
    lie2->enable_all_to_all_dump();
    lie2->set_output_dir(argv[3]);
    //lie2->set_output_dir(argv[4]);
    lie2->enable_IO();
    lie2->set_name(argv[0]);
    lie2->set_comm(mcomm);
    lie2->set_batch_size(10);
    lie2->execute();
    lie2->print_all_relation_size();
    delete lie2;
    }


    {
    relation* T1 = new relation(1, true, 2, 257, "T", "../data/g5955/path_2_1_2", FULL);
    relation* G1 = new relation(1, true, 2, 256, "G", argv[1], FULL);
    RAM* join1 = new RAM(true, 1);
    join1->add_relation(G1, false);
    join1->add_relation(T1, true);
    join1->add_rule(new parallel_join(T1, T1, DELTA, G1, FULL, {3, 1}));
    RAM* copy1 = new RAM(false, 0);
    copy1->add_relation(G1, false);
    copy1->add_relation(T1, true);
    copy1->add_rule(new parallel_copy(T1, G1, FULL, {1, 0}));
    LIE* lie1 = new LIE();
    lie1->add_relation(T1);
    lie1->add_relation(G1);
    lie1->add_scc(join1);
    lie1->add_scc(copy1);
    lie1->add_scc_dependance(copy1, join1);
    //lie1->set_sloav_mode(atoi(argv[3]));
    lie1->set_sloav_mode(0);
    lie1->enable_all_to_all_dump();
    lie1->set_output_dir(argv[2]);
    lie1->enable_IO();
    lie1->set_name(argv[0]);
    lie1->set_comm(mcomm);
    lie1->set_batch_size(10);
    lie1->execute();
    lie1->print_all_relation_size();
    delete lie1;
    }


    {
    relation* T2 = new relation(1, true, 2, 257, "T", "../data/g5955/path_2_1_2", FULL);
    relation* G2 = new relation(1, true, 2, 256, "G", argv[1], FULL);
    RAM* join2 = new RAM(true, 1);
    join2->add_relation(G2, false);
    join2->add_relation(T2, true);
    join2->add_rule(new parallel_join(T2, T2, DELTA, G2, FULL, {3, 1}));
    RAM* copy2 = new RAM(false, 0);
    copy2->add_relation(G2, false);
    copy2->add_relation(T2, true);
    copy2->add_rule(new parallel_copy(T2, G2, FULL, {1, 0}));
    LIE* lie2 = new LIE();
    lie2->add_relation(T2);
    lie2->add_relation(G2);
    lie2->add_scc(join2);
    lie2->add_scc(copy2);
    lie2->add_scc_dependance(copy2, join2);
    //lie2->set_sloav_mode(atoi(argv[5]));
    lie2->set_sloav_mode(1);
    lie2->enable_all_to_all_dump();
    lie2->set_output_dir(argv[3]);
    //lie2->set_output_dir(argv[4]);
    lie2->enable_IO();
    lie2->set_name(argv[0]);
    lie2->set_comm(mcomm);
    lie2->set_batch_size(10);
    lie2->execute();
    lie2->print_all_relation_size();
    delete lie2;
    }


    {
    relation* T1 = new relation(1, true, 2, 257, "T", "../data/g5955/path_2_1_2", FULL);
    relation* G1 = new relation(1, true, 2, 256, "G", argv[1], FULL);
    RAM* join1 = new RAM(true, 1);
    join1->add_relation(G1, false);
    join1->add_relation(T1, true);
    join1->add_rule(new parallel_join(T1, T1, DELTA, G1, FULL, {3, 1}));
    RAM* copy1 = new RAM(false, 0);
    copy1->add_relation(G1, false);
    copy1->add_relation(T1, true);
    copy1->add_rule(new parallel_copy(T1, G1, FULL, {1, 0}));
    LIE* lie1 = new LIE();
    lie1->add_relation(T1);
    lie1->add_relation(G1);
    lie1->add_scc(join1);
    lie1->add_scc(copy1);
    lie1->add_scc_dependance(copy1, join1);
    //lie1->set_sloav_mode(atoi(argv[3]));
    lie1->set_sloav_mode(0);
    lie1->enable_all_to_all_dump();
    lie1->set_output_dir(argv[2]);
    lie1->enable_IO();
    lie1->set_name(argv[0]);
    lie1->set_comm(mcomm);
    lie1->set_batch_size(10);
    lie1->execute();
    lie1->print_all_relation_size();
    delete lie1;
    }


    {
    relation* T2 = new relation(1, true, 2, 257, "T", "../data/g5955/path_2_1_2", FULL);
    relation* G2 = new relation(1, true, 2, 256, "G", argv[1], FULL);
    RAM* join2 = new RAM(true, 1);
    join2->add_relation(G2, false);
    join2->add_relation(T2, true);
    join2->add_rule(new parallel_join(T2, T2, DELTA, G2, FULL, {3, 1}));
    RAM* copy2 = new RAM(false, 0);
    copy2->add_relation(G2, false);
    copy2->add_relation(T2, true);
    copy2->add_rule(new parallel_copy(T2, G2, FULL, {1, 0}));
    LIE* lie2 = new LIE();
    lie2->add_relation(T2);
    lie2->add_relation(G2);
    lie2->add_scc(join2);
    lie2->add_scc(copy2);
    lie2->add_scc_dependance(copy2, join2);
    //lie2->set_sloav_mode(atoi(argv[5]));
    lie2->set_sloav_mode(1);
    lie2->enable_all_to_all_dump();
    lie2->set_output_dir(argv[3]);
    //lie2->set_output_dir(argv[4]);
    lie2->enable_IO();
    lie2->set_name(argv[0]);
    lie2->set_comm(mcomm);
    lie2->set_batch_size(10);
    lie2->execute();
    lie2->print_all_relation_size();
    delete lie2;
    }


    {
    relation* T1 = new relation(1, true, 2, 257, "T", "../data/g5955/path_2_1_2", FULL);
    relation* G1 = new relation(1, true, 2, 256, "G", argv[1], FULL);
    RAM* join1 = new RAM(true, 1);
    join1->add_relation(G1, false);
    join1->add_relation(T1, true);
    join1->add_rule(new parallel_join(T1, T1, DELTA, G1, FULL, {3, 1}));
    RAM* copy1 = new RAM(false, 0);
    copy1->add_relation(G1, false);
    copy1->add_relation(T1, true);
    copy1->add_rule(new parallel_copy(T1, G1, FULL, {1, 0}));
    LIE* lie1 = new LIE();
    lie1->add_relation(T1);
    lie1->add_relation(G1);
    lie1->add_scc(join1);
    lie1->add_scc(copy1);
    lie1->add_scc_dependance(copy1, join1);
    //lie1->set_sloav_mode(atoi(argv[3]));
    lie1->set_sloav_mode(0);
    lie1->enable_all_to_all_dump();
    lie1->set_output_dir(argv[2]);
    lie1->enable_IO();
    lie1->set_name(argv[0]);
    lie1->set_comm(mcomm);
    lie1->set_batch_size(10);
    lie1->execute();
    lie1->print_all_relation_size();
    delete lie1;
    }


    {
    relation* T2 = new relation(1, true, 2, 257, "T", "../data/g5955/path_2_1_2", FULL);
    relation* G2 = new relation(1, true, 2, 256, "G", argv[1], FULL);
    RAM* join2 = new RAM(true, 1);
    join2->add_relation(G2, false);
    join2->add_relation(T2, true);
    join2->add_rule(new parallel_join(T2, T2, DELTA, G2, FULL, {3, 1}));
    RAM* copy2 = new RAM(false, 0);
    copy2->add_relation(G2, false);
    copy2->add_relation(T2, true);
    copy2->add_rule(new parallel_copy(T2, G2, FULL, {1, 0}));
    LIE* lie2 = new LIE();
    lie2->add_relation(T2);
    lie2->add_relation(G2);
    lie2->add_scc(join2);
    lie2->add_scc(copy2);
    lie2->add_scc_dependance(copy2, join2);
    //lie2->set_sloav_mode(atoi(argv[5]));
    lie2->set_sloav_mode(1);
    lie2->enable_all_to_all_dump();
    lie2->set_output_dir(argv[3]);
    //lie2->set_output_dir(argv[4]);
    lie2->enable_IO();
    lie2->set_name(argv[0]);
    lie2->set_comm(mcomm);
    lie2->set_batch_size(10);
    lie2->execute();
    lie2->print_all_relation_size();
    delete lie2;
    }



    /*
    relation* T3 = new relation(1, true, 2, 257, "T", "../data/g5955/path_2_1_2", FULL);
    relation* G3 = new relation(1, true, 2, 256, "G", argv[1], FULL);
    RAM* join3 = new RAM(true, 1);
    join3->add_relation(G3, false);
    join3->add_relation(T3, true);
    join3->add_rule(new parallel_join(T3, T3, DELTA, G3, FULL, {3, 1}));
    RAM* copy3 = new RAM(false, 0);
    copy3->add_relation(G3, false);
    copy3->add_relation(T3, true);
    copy3->add_rule(new parallel_copy(T3, G3, FULL, {1, 0}));
    LIE* lie3 = new LIE();
    lie3->add_relation(T3);
    lie3->add_relation(G3);
    lie3->add_scc(join3);
    lie3->add_scc(copy3);
    lie3->add_scc_dependance(copy3, join3);
    lie3->set_sloav_mode(3);
    lie3->enable_all_to_all_dump();
    lie3->set_output_dir(argv[4]);
    lie3->enable_IO();
    lie3->set_name(argv[0]);
    lie3->set_comm(mcomm);
    lie3->set_batch_size(10);
    lie3->execute();
    lie3->print_all_relation_size();
    delete lie3;

    relation* T4 = new relation(1, true, 2, 257, "T", "../data/g5955/path_2_1_2", FULL);
    relation* G4 = new relation(1, true, 2, 256, "G", argv[1], FULL);
    RAM* join4 = new RAM(true, 1);
    join4->add_relation(G4, false);
    join4->add_relation(T4, true);
    join4->add_rule(new parallel_join(T4, T4, DELTA, G4, FULL, {3, 1}));
    RAM* copy4 = new RAM(false, 0);
    copy4->add_relation(G4, false);
    copy4->add_relation(T4, true);
    copy4->add_rule(new parallel_copy(T4, G4, FULL, {1, 0}));
    LIE* lie4 = new LIE();
    lie4->add_relation(T4);
    lie4->add_relation(G4);
    lie4->add_scc(join4);
    lie4->add_scc(copy4);
    lie4->add_scc_dependance(copy4, join4);
    lie4->set_sloav_mode(4);
    lie4->enable_all_to_all_dump();
    lie4->set_output_dir(argv[5]);
    lie4->enable_IO();
    lie4->set_name(argv[0]);
    lie4->set_comm(mcomm);
    lie4->set_batch_size(10);
    lie4->execute();
    lie4->print_all_relation_size();
    delete lie4;
    */

    mcomm.destroy();
    return 0;
}
