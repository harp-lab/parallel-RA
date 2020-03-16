#include "balanced_hash_relation.h"
#include "comm.h"
#include "parallel_join.h"
#include "RA_tasks.h"


int main(int argc, char **argv)
{
    mpi_comm mcomm;
    mcomm.create(argc, argv);

    relation* G0 = new relation();
    relation* G1 = new relation();
    relation* T0 = new relation();
    relation* T1 = new relation();
    relation* T4 = new relation();
    relation* W1 = new relation();
    relation* W2 = new relation();

    G0->initialize_relation(2, mcomm);
    G1->initialize_relation(2, mcomm);
    T0->initialize_relation(2, mcomm);
    T1->initialize_relation(2, mcomm);
    T4->initialize_relation(2, mcomm);
    W1->initialize_relation(2, mcomm);
    W2->initialize_relation(2, mcomm);


    G0->read_from_file(argv[1], FULL);


    if (mcomm.get_rank() == 0)
        std::cout << "-----------------------" << std::endl << std::endl;

    parallel_copy* copy_1 = new parallel_copy();
    copy_1->set_comm(mcomm);
    copy_1->set_copy_input(G0, FULL);
    copy_1->set_copy_output(G1);
    int rename_index_array_copy1[2] = {1, 0};
    int rename_index_array_size_copy1 = 2;
    copy_1->set_copy_rename_index(rename_index_array_copy1, rename_index_array_size_copy1);

    RAM scc0;
    scc0.push_relation(G0);
    scc0.push_relation(G1);
    scc0.set_comm(mcomm);
    scc0.push_back(copy_1);
    scc0.execute();

    if (mcomm.get_rank() == 0)
        std::cout << "-----------------------" << std::endl << std::endl;


    G1->read_from_relation(G1, DELTA);
    G1->flush_full();
    parallel_join* join_1 = new parallel_join();
    join_1->set_comm(mcomm);
    join_1->set_join_input0(G0, FULL);
    join_1->set_join_input1(G1, DELTA);
    join_1->set_join_output(W1);
    int project_and_rename_index_array_join1[3] = {-1, 0, 1};
    int project_and_rename_index_array_size_join1 = 3;
    join_1->set_join_projection_index(project_and_rename_index_array_join1, project_and_rename_index_array_size_join1);
    join_1->set_join_column_count(1);


    RAM scc1;
    scc1.push_relation(G0);
    scc1.push_relation(G1);
    scc1.push_relation(W1);
    scc1.set_comm(mcomm);
    scc1.push_back(join_1);
    scc1.execute();

    if (mcomm.get_rank() == 0)
        std::cout << "-----------------------" << std::endl << std::endl;


    parallel_copy* copy_2 = new parallel_copy();
    copy_2->set_comm(mcomm);
    copy_2->set_copy_input(G0, FULL);
    copy_2->set_copy_output(T0);
    int rename_index_array_copy2[2] = {0, 1};
    int rename_index_array_size_copy2 = 2;
    copy_2->set_copy_rename_index(rename_index_array_copy2, rename_index_array_size_copy2);

    RAM scc2;
    scc2.push_relation(G0);
    scc2.push_relation(T0);
    scc2.set_comm(mcomm);
    scc2.push_back(copy_2);
    scc2.execute();

    if (mcomm.get_rank() == 0)
        std::cout << "-----------------------" << std::endl << std::endl;

    parallel_copy* copy_3 = new parallel_copy();
    copy_3->set_comm(mcomm);
    copy_3->set_copy_input(W1, FULL);
    copy_3->set_copy_output(T0);
    int rename_index_array_copy3[2] = {1, 0};
    int rename_index_array_size_copy3 = 2;
    copy_3->set_copy_rename_index(rename_index_array_copy3, rename_index_array_size_copy3);

    RAM scc3;
    scc3.push_relation(W1);
    scc3.push_relation(T0);
    scc3.set_comm(mcomm);
    scc3.push_back(copy_3);
    scc3.execute();

    if (mcomm.get_rank() == 0)
        std::cout << "-----------------------" << std::endl << std::endl;

    parallel_copy* copy_4 = new parallel_copy();
    copy_4->set_comm(mcomm);
    copy_4->set_copy_input(T0, FULL);
    copy_4->set_copy_output(T1);
    int rename_index_array_copy4[2] = {1, 0};
    int rename_index_array_size_copy4 = 2;
    copy_4->set_copy_rename_index(rename_index_array_copy4, rename_index_array_size_copy4);

    RAM scc4;
    scc4.push_relation(T0);
    scc4.push_relation(T1);
    scc4.set_comm(mcomm);
    scc4.push_back(copy_4);
    scc4.execute();

    if (mcomm.get_rank() == 0)
        std::cout << "-----------------------" << std::endl << std::endl;

    T1->read_from_relation(T1, DELTA);
    T1->flush_full();
    parallel_join* join_2 = new parallel_join();
    join_2->set_comm(mcomm);
    join_2->set_join_input0(T0, FULL);
    join_2->set_join_input1(T1, DELTA);
    join_2->set_join_output(W2);
    int project_and_rename_index_array_join2[3] = {-1, 0, 1};
    int project_and_rename_index_array_size_join2 = 3;
    join_2->set_join_projection_index(project_and_rename_index_array_join2, project_and_rename_index_array_size_join2);
    join_2->set_join_column_count(1);


    RAM scc5;
    scc5.push_relation(T0);
    scc5.push_relation(T1);
    scc5.push_relation(W2);
    scc5.set_comm(mcomm);
    scc5.push_back(join_2);
    scc5.execute();

    if (mcomm.get_rank() == 0)
        std::cout << "-----------------------" << std::endl << std::endl;


    parallel_copy* copy_5 = new parallel_copy();
    copy_5->set_comm(mcomm);
    copy_5->set_copy_input(G0, FULL);
    copy_5->set_copy_output(T4);
    int rename_index_array_copy5[2] = {0, 1};
    int rename_index_array_size_copy5 = 2;
    copy_5->set_copy_rename_index(rename_index_array_copy5, rename_index_array_size_copy5);

    RAM scc6;
    scc6.push_relation(G0);
    scc6.push_relation(T4);
    scc6.set_comm(mcomm);
    scc6.push_back(copy_5);
    scc6.execute();

    if (mcomm.get_rank() == 0)
        std::cout << "-----------------------" << std::endl << std::endl;

    parallel_copy* copy_6 = new parallel_copy();
    copy_6->set_comm(mcomm);
    copy_6->set_copy_input(W2, FULL);
    copy_6->set_copy_output(T4);
    int rename_index_array_copy6[2] = {1, 0};
    int rename_index_array_size_copy6 = 2;
    copy_6->set_copy_rename_index(rename_index_array_copy6, rename_index_array_size_copy6);

    RAM scc7;
    scc7.push_relation(W2);
    scc7.push_relation(T4);
    scc7.set_comm(mcomm);
    scc7.push_back(copy_6);
    scc7.execute();

    if (mcomm.get_rank() == 0)
        std::cout << "-----------------------" << std::endl << std::endl;


    T4->read_from_relation(T4, DELTA);
    T4->flush_full();
    parallel_join* join_3 = new parallel_join();
    join_3->set_comm(mcomm);
    join_3->set_join_input0(T4, DELTA);
    join_3->set_join_input1(W2, FULL);
    join_3->set_join_output(T4);
    int project_and_rename_index_array_join3[3] = {-1, 0, 1};
    int project_and_rename_index_array_size_join3 = 3;
    join_3->set_join_projection_index(project_and_rename_index_array_join3, project_and_rename_index_array_size_join3);
    join_3->set_join_column_count(1);


    RAM scc8;
    scc8.push_relation(T4);
    scc8.push_relation(W2);
    scc8.set_comm(mcomm);
    scc8.push_back(join_3);
    scc8.execute();

    if (mcomm.get_rank() == 0)
        std::cout << "-----------------------" << std::endl << std::endl;


    delete copy_1;
    delete copy_2;
    delete copy_3;
    delete copy_4;
    delete copy_5;
    delete copy_6;
    delete join_1;
    delete join_2;
    delete join_3;


    delete G0;
    delete G1;
    delete W1;
    delete W2;
    delete T0;
    delete T1;
    delete T4;

    mcomm.destroy();

    return 0;
}
