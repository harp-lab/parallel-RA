#include "../src/parallel_RA_inc.h"


int main(int argc, char **argv)
{
    mpi_comm mcomm;
    mcomm.create(argc, argv);

    relation* rel_var_ref_1_1 = new relation(1, 1, 258, "/var/tmp/g5394/var-ref_1_1", FULL);
    relation* rel_app_2_1_2 = new relation(2, 2, 256, "/var/tmp/g5394/app_2_1_2", FULL);
    relation* rel_inter_body24_3_2 = new relation(1, 3, 257, "/var/tmp/g5394/inter-body24_3_2", FULL);
    relation* rel_inter_body24_3_1_2_3 = new relation(3, 3, 257, "/var/tmp/g5394/inter-body24_3_1_2_3", FULL);
    relation* rel_args_1_1 = new relation(1, 1, 261, "/var/tmp/g5394/args_1_1", FULL);
    relation* rel_app_2_1 = new relation(1, 2, 256, "/var/tmp/g5394/app_2_1", FULL);
    relation* rel_app_2_2 = new relation(1, 2, 256, "/var/tmp/g5394/app_2_2", FULL);
    relation* rel_lambda_2_2 = new relation(1, 2, 260, "/var/tmp/g5394/lambda_2_2", FULL);
    relation* rel_lambda_2_1_2 = new relation(2, 2, 260, "/var/tmp/g5394/lambda_2_1_2", FULL);
    relation* rel_free_2_1_2 = new relation(2, 2, 259, "/var/tmp/g5394/free_2_1_2", FULL);
    relation* rel_args_1_0 = new relation(1, 1, 261, "/var/tmp/g5394/args_1_0", FULL);
    relation* rel_free_2_2 = new relation(1, 2, 259, "/var/tmp/g5394/free_2_2", FULL);
    relation* rel_inter_body26_3_1_2_3 = new relation(3, 3, 262, "/var/tmp/g5394/inter-body26_3_1_2_3", FULL);
    relation* rel_var_ref_1_ = new relation(0, 1, 258, "/var/tmp/g5394/var-ref_1_", FULL);
    relation* rel_inter_body26_3_3_2 = new relation(2, 3, 262, "/var/tmp/g5394/inter-body26_3_3_2", FULL);

    RAM* scc5395 = new RAM(false);
    scc5395->add_relation(rel_app_2_1, true);
    scc5395->add_relation(rel_app_2_1_2, true);
    scc5395->add_rule(new parallel_acopy(rel_app_2_1, rel_app_2_1_2, DELTA, {0, 2, 1}));

    RAM* scc5396 = new RAM(false);
    scc5396->add_relation(rel_app_2_2, true);
    scc5396->add_relation(rel_app_2_1_2, true);
    scc5396->add_rule(new parallel_acopy(rel_app_2_2, rel_app_2_1_2, DELTA, {2, 0, 1}));

    RAM* scc5397 = new RAM(false);
    scc5397->add_relation(rel_free_2_1_2, true);
    scc5397->add_relation(rel_var_ref_1_1, false);
    scc5397->add_rule(new parallel_copy(rel_free_2_1_2, rel_var_ref_1_1, FULL, {0, 1}));

    RAM* scc5398 = new RAM(false);
    scc5398->add_relation(rel_lambda_2_1_2, true);
    scc5398->add_relation(rel_lambda_2_2, true);
    scc5398->add_rule(new parallel_acopy(rel_lambda_2_2, rel_lambda_2_1_2, DELTA, {2, 0, 1}));

    RAM* scc5399 = new RAM(true);
    scc5399->add_relation(rel_inter_body26_3_3_2, true);
    scc5399->add_relation(rel_inter_body26_3_1_2_3, true);
    scc5399->add_relation(rel_free_2_2, true);
    scc5399->add_relation(rel_args_1_0, false);
    scc5399->add_relation(rel_free_2_1_2, true);
    scc5399->add_relation(rel_lambda_2_2, false);
    scc5399->add_relation(rel_app_2_2, false);
    scc5399->add_relation(rel_app_2_1, false);
    scc5399->add_relation(rel_inter_body24_3_1_2_3, true);
    scc5399->add_relation(rel_inter_body24_3_2, true);
    scc5399->add_rule(new parallel_acopy(rel_free_2_2, rel_free_2_1_2, DELTA, {2, 0, 1}));
    scc5399->add_rule(new parallel_join(rel_free_2_1_2, rel_app_2_2, FULL, rel_free_2_2, DELTA, {-1, 1, -1, -1, 0}));
    scc5399->add_rule(new parallel_join(rel_free_2_1_2, rel_app_2_1, FULL, rel_free_2_2, DELTA, {-1, 1, -1, -1, 0}));
    scc5399->add_rule(new parallel_acopy(rel_inter_body24_3_2, rel_inter_body24_3_1_2_3, DELTA, {2, 0, 3, 1}));
    scc5399->add_rule(new parallel_acopy(rel_inter_body26_3_3_2, rel_inter_body26_3_1_2_3, DELTA, {3, 0, 1, 2}));
    scc5399->add_rule(new parallel_copy_filter(rel_free_2_1_2, rel_inter_body26_3_3_2, DELTA, {0, -1, -1, 1}, [](const u64* const data){ return !(data[0] == data[1]); }));
    scc5399->add_rule(new parallel_join(rel_inter_body24_3_1_2_3, rel_free_2_2, DELTA, rel_lambda_2_2, FULL, {-1, -1, 2, 0, 1}));
    scc5399->add_rule(new parallel_join(rel_inter_body26_3_1_2_3, rel_args_1_0, FULL, rel_inter_body24_3_2, DELTA, {-1, 1, -1, 0, 2}));

    RAM* scc5400 = new RAM(false);
    scc5400->add_relation(rel_args_1_0, true);
    scc5400->add_relation(rel_args_1_1, true);
    scc5400->add_rule(new parallel_acopy(rel_args_1_0, rel_args_1_1, DELTA, {1, 0}));

    RAM* scc5401 = new RAM(false);
    scc5401->add_relation(rel_var_ref_1_, true);
    scc5401->add_relation(rel_var_ref_1_1, true);
    scc5401->add_rule(new parallel_acopy(rel_var_ref_1_, rel_var_ref_1_1, DELTA, {1, 0}));

    LIE* lie = new LIE();
    lie->add_relation(rel_var_ref_1_1);
    lie->add_relation(rel_app_2_1_2);
    lie->add_relation(rel_inter_body24_3_2);
    lie->add_relation(rel_inter_body24_3_1_2_3);
    lie->add_relation(rel_args_1_1);
    lie->add_relation(rel_app_2_1);
    lie->add_relation(rel_app_2_2);
    lie->add_relation(rel_lambda_2_2);
    lie->add_relation(rel_lambda_2_1_2);
    lie->add_relation(rel_free_2_1_2);
    lie->add_relation(rel_args_1_0);
    lie->add_relation(rel_free_2_2);
    lie->add_relation(rel_inter_body26_3_1_2_3);
    lie->add_relation(rel_var_ref_1_);
    lie->add_relation(rel_inter_body26_3_3_2);
    lie->add_scc(scc5395);
    lie->add_scc(scc5396);
    lie->add_scc(scc5397);
    lie->add_scc(scc5398);
    //lie->add_scc(scc5399);
    lie->add_scc(scc5400);
    lie->add_scc(scc5401);
    //lie->add_scc_dependance(scc5395, scc5399);
    //lie->add_scc_dependance(scc5396, scc5399);
    //lie->add_scc_dependance(scc5397, scc5399);
    //lie->add_scc_dependance(scc5398, scc5399);
    //lie->add_scc_dependance(scc5400, scc5399);



    lie->set_comm(mcomm);
    lie->set_batch_size(1);


    lie->execute();

    rel_free_2_2->print();
    rel_app_2_1->print();
    rel_app_2_2->print();
    rel_app_2_1_2->print();


    mcomm.destroy();
    return 0;
}
