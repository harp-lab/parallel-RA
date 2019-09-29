#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <mpi.h>

#include "btree.h"
#include "btree_relation.h"

#include "comm.h"
#include "parallel_io.h"
#include "RA.h"


int main(int argc, char **argv)
{
    mpi_comm mcomm;
    mcomm.create(argc, argv);

    parallel_io ioG(argv[1]);
    ioG.parallel_read_input_relation_from_file_to_local_buffer(mcomm.get_rank(), mcomm.get_nprocs());
    ioG.buffer_data_to_hash_buffer(0 /*hash_column_index*/, mcomm.get_nprocs(), mcomm.get_comm());
    ioG.delete_raw_buffers();

    parallel_io ioT(argv[1]);
    ioT.parallel_read_input_relation_from_file_to_local_buffer(mcomm.get_rank(), mcomm.get_nprocs());
    ioT.buffer_data_to_hash_buffer(1 /*hash_column_index*/, mcomm.get_nprocs(), mcomm.get_comm());
    ioT.delete_raw_buffers();

    relation<2> T(ioT.get_col_count() * ioT.get_row_count(), ioT.get_hash_buffer());
    relation<2> G(ioG.get_col_count() * ioG.get_row_count(), ioG.get_hash_buffer());

    relation<2> * dT = new relation<2>;
    dT->initialize(ioT.get_col_count() * ioT.get_row_count(), ioT.get_hash_buffer());

    ioT.delete_hash_buffers();
    ioG.delete_hash_buffers();

    RA rj;
    int lb = 0;
    int running_t_count = ioT.get_row_count();
    double time = 0;

    dT = rj.parallel_join(dT, G, T, 0, &lb, &running_t_count, &time, mcomm.get_nprocs(), mcomm.get_rank(), mcomm.get_comm());

    int lc = 1;
    while(true)
    {
      dT = rj.parallel_join(dT, G, T, lc, &lb, &running_t_count, &time, mcomm.get_nprocs(), mcomm.get_rank(), mcomm.get_comm());

      if (lb == 1)  break;
      lc++;
    }
    delete dT;


    mcomm.destroy();

    return 0;
}







