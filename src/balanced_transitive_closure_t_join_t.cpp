#include "balanced_hash_relation.h"
#include "comm.h"
#include "parallel_join.h"
#include "RA_tasks.h"


int main(int argc, char **argv)
{

    mpi_comm mcomm;
    mcomm.create(argc, argv);
    mcomm.set_number_of_buckets(1);

    u32 sub_buckets_per_bucket_T = 1;

    u32 arity = 2;
    u32 join_column_count = 1;


    relation* T0 = new relation();
    int rename_and_project_copy0[2] = {0,1};
    T0->initialize_with_rename_and_projection(arity, join_column_count, DYNAMIC, sub_buckets_per_bucket_T, argv[1], DELTA, rename_and_project_copy0, mcomm);

    relation* T1 = new relation();
    int rename_and_project_copy1[2] = {1,0};
    T1->initialize_with_rename_and_projection(arity, join_column_count, DYNAMIC, sub_buckets_per_bucket_T, argv[1], DELTA, rename_and_project_copy1, mcomm);


#if 0
    /////// SCC 0
    int rename_and_project_copy0[2] = {0,1};
    parallel_copy copy_0;
    copy_0.input(G, FULL);
    copy_0.output(T0);
    copy_0.projection_index(rename_and_project_copy0);

    int rename_and_project_copy1[2] = {1,0};
    parallel_copy copy_1;
    copy_1.input(G, FULL);
    copy_1.output(T1);
    copy_1.projection_index(rename_and_project_copy1);

    RAM scc0;
    scc0.push_back(copy_0);
    scc0.push_back(copy_1);
    scc0.execute();
#endif


    /////// SCC 1
    int rename_and_project_join[3] = {-1,0,1};
    parallel_RA join_0(JOIN, mcomm);
    join_0.join_input0(T1, DELTA);
    join_0.join_input1(T0, DELTA);
    join_0.join_output(T0);
    join_0.set_projection_index(rename_and_project_join);

    parallel_RA join_1(JOIN, mcomm);
    join_1.join_input0(T1, FULL);
    join_1.join_input1(T0, DELTA);
    join_1.join_output(T0);
    join_1.set_projection_index(rename_and_project_join);

    parallel_RA join_2(JOIN, mcomm);
    join_2.join_input0(T1, DELTA);
    join_2.join_input1(T0, FULL);
    join_2.join_output(T0);
    join_2.set_projection_index(rename_and_project_join);

    int rename_and_project_acopy[2] = {1,0};
    parallel_RA copy_0(COPY, mcomm);
    copy_0.copy_input(T0, DELTA);
    copy_0.copy_output(T1);
    copy_0.set_projection_index(rename_and_project_acopy);

    RAM scc1(mcomm);
    scc1.push_back(join_0);
    scc1.push_back(join_1);
    scc1.push_back(join_2);
    scc1.push_back(copy_0);
    scc1.execute();

    delete T0;
    delete T1;

    mcomm.destroy();

    return 0;
}
