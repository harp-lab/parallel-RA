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

    relation* G0 = new relation();
    relation* G1 = new relation();
    relation* W = new relation();
    relation* T = new relation();

    relation* Tnew = new relation();

    double t1 = MPI_Wtime();
    int rename_and_project_copy[2] = {1,0};
    G0->initialize(arity, join_column_count, STATIC, sub_buckets_per_bucket_G, argv[1], FULL, mcomm, 0);
    W->initialize_empty(2, join_column_count, sub_buckets_per_bucket_G, mcomm);
    T->initialize_empty(2, join_column_count, sub_buckets_per_bucket_G, mcomm);
    G1->initialize_with_rename_and_projection(arity, join_column_count, DYNAMIC, sub_buckets_per_bucket_T, argv[1], FULL, rename_and_project_copy, mcomm, 1);
    double t2 = MPI_Wtime();

    if (mcomm.get_rank() == 0)
        std::cout << "----------------- [1] Time " << t2 - t1 << std::endl;

    //G0->print();
    //G1->print();


    double t3 = MPI_Wtime();
    parallel_RA* join_1 = new parallel_RA(JOIN, mcomm);
    join_1->set_RA_typex(0);
    join_1->join_input0(G0, FULL);
    join_1->join_input1(G1, FULL);
    join_1->join_output(W);
    join_1->set_join_projection_index(-1, 0, 1);

    RAM scc1(mcomm);
    scc1.push_back(join_1);
    scc1.set_threshold(atoi(argv[2]));
    scc1.set_refinement_interval(atoi(argv[3]));
    scc1.set_refinement_factor(atof(argv[4]));
    scc1.set_refinement_chooser(atoi(argv[5]));
    scc1.execute();
    double t4 = MPI_Wtime();

    if (mcomm.get_rank() == 0)
        std::cout << "----------------- [2] Time " << t4 - t3 << std::endl;

    //W->print();



    double t5 = MPI_Wtime();
    /////// SCC 2
    parallel_RA* copy_1 = new parallel_RA(COPY, mcomm);
    copy_1->copy_input(G0, FULL);
    copy_1->copy_output(T);
    copy_1->set_copy_projection_index(0, 1);
    RAM scc2(mcomm);
    scc2.push_back(copy_1);
    scc2.set_threshold(atoi(argv[2]));
    scc2.set_refinement_interval(atoi(argv[3]));
    scc2.set_refinement_factor(atof(argv[4]));
    scc2.set_refinement_chooser(atoi(argv[5]));
    scc2.execute();
    double t6 = MPI_Wtime();

    if (mcomm.get_rank() == 0)
        std::cout << "----------------- [3] Time " << t6 - t5 << std::endl;

    //T->print();



    double t7 = MPI_Wtime();
    /////// SCC 3
    parallel_RA* copy_2 = new parallel_RA(COPY, mcomm);
    copy_2->copy_input(W, FULL);
    copy_2->copy_output(T);
    copy_2->set_copy_projection_index(1, 0);
    RAM scc3(mcomm);
    scc3.push_back(copy_2);
    scc3.set_threshold(atoi(argv[2]));
    scc3.set_refinement_interval(atoi(argv[3]));
    scc3.set_refinement_factor(atof(argv[4]));
    scc3.set_refinement_chooser(atoi(argv[5]));
    scc3.execute();
    double t8 = MPI_Wtime();

    Tnew->initialize_delta_with_relation(arity, join_column_count, sub_buckets_per_bucket_G, mcomm, T);

    //Tnew->print();
    //W->print();

    if (mcomm.get_rank() == 0)
        std::cout << "----------------- [4] Time " << t8 - t7 << std::endl;


    double t9 = MPI_Wtime();
    /////// SCC 4
    parallel_RA* join_2 = new parallel_RA(JOIN, mcomm);
    join_2->set_RA_typex(1);
    join_2->join_input0(Tnew, DELTA);
    join_2->join_output(Tnew);
    join_2->join_input1(W, FULL);
    join_2->set_join_projection_index(-1, 0, 1);
    RAM scc4(mcomm);
    scc4.push_back(join_2);
    scc4.set_threshold(atoi(argv[2]));
    scc4.set_refinement_interval(atoi(argv[3]));
    scc4.set_refinement_factor(atof(argv[4]));
    scc4.set_refinement_chooser(atoi(argv[5]));
    scc4.execute();
    double t10 = MPI_Wtime();

    if (mcomm.get_rank() == 0)
        std::cout << "----------------- [5] Time " << t10 - t9 << std::endl;

    //Tnew->print();

    delete copy_1;
    delete copy_2;
    delete join_1;
    delete join_2;


    delete Tnew;
    delete G0;
    delete G1;
    delete W;
    delete T;

    mcomm.destroy();

    return 0;
}
