#include "parallel_RA_inc.h"



int main(int argc, char **argv)
{
    mpi_comm mcomm;
    mcomm.create(argc, argv);

    u32 batch_size = atoi(argv[1]);
    double batch_time = atof(argv[2]);
    float threshold = atof(argv[3]);
    int mode = atoi(argv[4]);
    u32 task_count = atoi(argv[5]);

    LIE engine;
    engine.set_comm(mcomm);
    engine.set_task_threshold(threshold);

    engine.set_batch_time(batch_time);
    engine.set_batch_size(batch_size);

    engine.set_mode(mode);
    std::vector<std::vector<RAM>> taskgraph(2, std::vector<RAM>(task_count));

    relation* G = new relation[task_count];
    relation* T = new relation[task_count];

    RAM* scc0 = new RAM[task_count];
    RAM* scc1 = new RAM[task_count];

    for (u32 i=0; i< task_count; i++)
    {
        G[i].initialize(2, atoi(argv[6+2*i]), argv[7+2*i], FULL);
        T[i].initialize(2);

        scc0[i].push_relation(&G[i]);
        scc0[i].push_relation(&T[i]);
        scc0[i].set_iteration_count(1);
        scc0[i].push_rule(new parallel_copy(&G[i], FULL, &T[i], {1, 0}));
        taskgraph[0][i] = scc0[i];

        scc1[i].push_relation(&(G[i]));
        scc1[i].push_relation(&(T[i]));
        scc1[i].allow_print_result();
        scc1[i].push_rule(new parallel_join(&G[i], FULL, &T[i], DELTA, &T[i], 1, {-1, 0, 1}));
        taskgraph[1][i] = scc1[i];
    }
    engine.set_taskgraph(taskgraph);
    engine.execute();

    delete[] scc0;
    delete[] scc1;
    delete[] G;
    delete[] T;

    mcomm.destroy();
    return 0;
}
