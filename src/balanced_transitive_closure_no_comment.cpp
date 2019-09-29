#include "balanced_hash_relation.h"
#include "comm.h"
#include "balanced_RA.h"

int main(int argc, char **argv)
{
    mpi_comm mcomm;
    mcomm.create(argc, argv);

    u32 buckets = mcomm.get_nprocs();
    u32 sub_buckets_per_bucket = 1;

    relation<2> G;
    G.initialize(buckets, sub_buckets_per_bucket, argv[1], FULL, mcomm);

    relation<2> T;
    T.initialize(buckets, sub_buckets_per_bucket, argv[1], FULL_AND_DELTA, mcomm);

    double lb_factor = atof(4);
    G.load_balance(lb_factor);
    T.load_balance(lb_factor);

    balanced_RA<2,2,2> b_ra;

    int project_indices[3] = {0,1,1};
    int renaming_indices[2] = {0,1};
    
    while(b_ra.parallel_map_join(G, FULL, T, DELTA, 1, T, FULL_AND_DELTA, project_indices, renaming_indices))
        T.load_balance(lb_factor);

    mcomm.destroy();

    return 0;
}
