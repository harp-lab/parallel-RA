#include <chrono>
#include <sys/stat.h>
#include <errno.h>
#include <limits.h>
#include <arpa/inet.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <unistd.h>
#include <iostream>
#include <fstream>

#include "btree.h"
#include "relation.h"
#include "RA.h"


static void read_input_relation_from_file_to_local_buffer(const char *file_name, u64** read_buffer, u32* row_count, u32* col_count)
{
    u32 r_count;
    u32 c_count;

    char meta_data_filename[1024];
    sprintf(meta_data_filename, "%s/meta_data.txt", file_name);

    FILE *fp_in;
    fp_in = fopen(meta_data_filename, "r");
    if (fscanf (fp_in, "(row count)\n%d\n(col count)\n%d", &r_count, &c_count) != 2)
    {
        printf("Wrong input format (Meta Data) %d %d\n", r_count, c_count);
        exit(0);
    }
    fclose(fp_in);

    *row_count = r_count;
    *col_count = c_count;


    char data_filename[1024];
    sprintf(data_filename, "%s/data.raw", file_name);
    int fp = open(data_filename, O_RDONLY);

    *read_buffer = new u64[r_count * c_count];
    u32 rb_size = (u32)read(fp, *read_buffer, r_count * c_count * sizeof(u64));
    if (rb_size != r_count * c_count * sizeof(u64))
    {
        printf("Wrong input format (DATA) %d %d\n", rb_size, (int)r_count * c_count * (int) sizeof(u64) );
        exit(0);
    }
    close(fp);

    return;
}


int main(int argc, char **argv)
{
    u64 *input_buffer = NULL;
    u32 entry_count;
    u32 col_count;
    read_input_relation_from_file_to_local_buffer(argv[1], &input_buffer, &entry_count, &col_count);

    relation<2> G(col_count * entry_count, input_buffer);
    relation<2> T;

    tuple<2> t;
    t[0] = -1; t[1] = -1;
    tuple<2> selectall(t);

    for (relation<2>::iter it(G, selectall); it.more(); it.advance())
    {
      tuple<2> t1;
      t1[0] = (*it)[0];
      t1[1] = (*it)[1];
      T.insert(t1);
    }

    int lc = 0;
    relation<2> delT = join(T, G, 0);
    union1(T, delT, 0);
    std::cout << "Loop count " << lc << std::endl;

    lc++;
    while(true)
    {
      std::cout << "Loop count " << lc << std::endl;
      delT = join(delT, G, lc);
      if (union1(T, delT, lc) == false)
        break;
      lc++;
    }

    char TCname[1024];
    sprintf(TCname, "%s_TC8", argv[1]);

    std::ofstream myfile;
    myfile.open (TCname);
    for (relation<2>::iter Tit(T, selectall); Tit.more(); Tit.advance())
      myfile << (*Tit)[0] << "\t" << (*Tit)[1] << "\n";
    myfile.close();

}
