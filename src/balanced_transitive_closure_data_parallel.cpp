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
    relation* T = new relation();
    relation* W = new relation();

    G0->initialize_relation(2, mcomm);
    G1->initialize_relation(2, mcomm);
    T->initialize_relation(2, mcomm);
    W->initialize_relation(2, mcomm);

    G0->read_from_file(argv[1], FULL);



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
    join_1->set_join_output(W);
    int project_and_rename_index_array_join1[3] = {-1, 0, 1};
    int project_and_rename_index_array_size_join1 = 3;
    join_1->set_join_projection_index(project_and_rename_index_array_join1, project_and_rename_index_array_size_join1);
    join_1->set_join_column_count(1);


    RAM scc1;
    scc1.push_relation(G0);
    scc1.push_relation(G1);
    scc1.push_relation(W);
    scc1.set_comm(mcomm);
    scc1.push_back(join_1);
    scc1.execute();

    if (mcomm.get_rank() == 0)
        std::cout << "-----------------------" << std::endl << std::endl;


    parallel_copy* copy_2 = new parallel_copy();
    copy_2->set_comm(mcomm);
    copy_2->set_copy_input(G0, FULL);
    copy_2->set_copy_output(T);
    int rename_index_array_copy2[2] = {0, 1};
    int rename_index_array_size_copy2 = 2;
    copy_2->set_copy_rename_index(rename_index_array_copy2, rename_index_array_size_copy2);

    RAM scc2;
    scc2.push_relation(G0);
    scc2.push_relation(T);
    scc2.set_comm(mcomm);
    scc2.push_back(copy_2);
    scc2.execute();

    if (mcomm.get_rank() == 0)
        std::cout << "-----------------------" << std::endl << std::endl;

    parallel_copy* copy_3 = new parallel_copy();
    copy_3->set_comm(mcomm);
    copy_3->set_copy_input(W, FULL);
    copy_3->set_copy_output(T);
    int rename_index_array_copy3[2] = {1, 0};
    int rename_index_array_size_copy3 = 2;
    copy_3->set_copy_rename_index(rename_index_array_copy3, rename_index_array_size_copy3);

    RAM scc3;
    scc3.push_relation(W);
    scc3.push_relation(T);
    scc3.set_comm(mcomm);
    scc3.push_back(copy_3);
    scc3.execute();

    if (mcomm.get_rank() == 0)
        std::cout << "-----------------------" << std::endl << std::endl;


    T->read_from_relation(T, DELTA);
    T->flush_full();
    parallel_join* join_2 = new parallel_join();
    join_2->set_comm(mcomm);
    join_2->set_join_input0(T, DELTA);
    join_2->set_join_input1(W, FULL);
    join_2->set_join_output(T);
    int project_and_rename_index_array_join2[3] = {-1, 0, 1};
    int project_and_rename_index_array_size_join2 = 3;
    join_2->set_join_projection_index(project_and_rename_index_array_join2, project_and_rename_index_array_size_join2);
    join_2->set_join_column_count(1);


    RAM scc4;
    scc4.push_relation(T);
    scc4.push_relation(W);
    scc4.set_comm(mcomm);
    scc4.push_back(join_2);
    //scc4.enable_logging();
    scc4.execute();

    //scc4.print_all_relation();

    if (mcomm.get_rank() == 0)
        std::cout << "-----------------------" << std::endl << std::endl;


    delete copy_1;
    delete copy_2;
    delete copy_3;
    delete join_1;
    delete join_2;


    delete G0;
    delete G1;
    delete W;
    delete T;

    mcomm.destroy();

    return 0;
}
