#include "balanced_hash_relation.h"
#include "comm.h"
#include "parallel_join.h"
#include "RA_tasks.h"


int main(int argc, char **argv)
{
    mpi_comm mcomm;
    mcomm.create(argc, argv);
    mcomm.set_number_of_buckets(1);

    int arity = 2;
    u32 sub_buckets_per_bucket_G = 1;
    u32 sub_buckets_per_bucket_T = 1;

    u32 join_column_count = 1;

    relation* G = new relation();
    relation* T = new relation();
    G->initialize(arity, join_column_count, STATIC, sub_buckets_per_bucket_G, argv[1], FULL, mcomm, 0);

    int rename_and_project_copy[2] = {1,0};
    T->initialize_with_rename_and_projection(arity, join_column_count, DYNAMIC, sub_buckets_per_bucket_T, argv[1], DELTA, rename_and_project_copy, mcomm, 1);


    /*
    parallel_RA copy_1(COPY, mcomm);
    copy_1.copy_input(G, FULL);
    copy_1.copy_output(T);
    copy_1.set_projection_index(rename_and_project_copy);

    RAM scc0(mcomm);
    scc0.push_back(copy_1);
    scc0.execute();
    */

    /////// SCC 1 (exchanged)
    //int rename_and_project_join[3] = {-1,0,1};
    parallel_RA* join_1 = new parallel_RA(JOIN, mcomm);
    join_1->join_input0(G, FULL);
    join_1->join_input1(T, DELTA);
    join_1->join_output(T);

    join_1->set_join_projection_index(-1, 0, 1);

#if 1
    RAM scc1(mcomm);
    scc1.push_back(join_1);
    scc1.set_threshold(atoi(argv[2]));
    scc1.set_refinement_interval(atoi(argv[3]));
    scc1.set_refinement_factor(atof(argv[4]));
    scc1.set_refinement_chooser(atoi(argv[5]));
    scc1.execute();
#endif

    delete G;
    delete T;
    delete join_1;

    mcomm.destroy();

    return 0;
}
