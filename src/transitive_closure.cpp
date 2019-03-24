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
#include "btree_relation.h"
//#include "RA.h"


uint64_t utime()
{
    return std::chrono::high_resolution_clock::now().time_since_epoch().count();
}

relation<2> * join(relation<2>* delT, relation<2>& G, relation<2>& T, int lc, int* lb, int* running_t_count, u64* running_time)
{
    u64 start = utime();
    tuple<2> t;
    t[0] = -1;
    t[1] = -1;
    tuple<2> selectall(t);
    relation<2> * delTT = new relation<2>;

    int count = 0;
    int tcount = 0;
    int tuple_count = 0;

    for (relation<2>::iter dit(*delT, selectall); dit.more(); dit.advance())
    {
      tuple<2> s;
      s[0] = (*dit)[1];
      s[1] = -1;
      tuple<2> select(s);

      for (relation<2>::iter git(G, select); git.more(); git.advance())
      {
        tuple<2> dt;
        dt[0] = (*dit)[0];
        dt[1] = (*git)[1];
        tuple_count++;

        if (T.insert(dt) == true)
        {
            tcount++;
            if (delTT->insert(dt) == true)
              count++;
        }
      }
    }

    if (count == 0)
        *lb = 1;
    else
        *lb = 0;

    *running_t_count = *running_t_count + tcount;

    u64 end = utime();
    u64 dTime = (end - start) / 1000000;
    *running_time = *running_time + dTime;

    std::cout << lc << " [" << dTime << "]  [" << *running_time << "] Tuple count " << tuple_count << " Delta count " << count << " T count: " << *running_t_count << " : " << std::endl;

    delete delT;

    return delTT;
}


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
    relation<2> T;// = new relation<2>;
    relation<2> * dT = new relation<2>;

    delete[] input_buffer;

    tuple<2> t;
    t[0] = -1; t[1] = -1;
    tuple<2> selectall(t);

    int running_t_count = 0;
    for (relation<2>::iter it(G, selectall); it.more(); it.advance())
    {
      tuple<2> t1;
      t1[0] = (*it)[0];
      t1[1] = (*it)[1];

      if (T.insert(t1) == true)
          running_t_count++;
      dT->insert(t1);
    }
    std::cout << "Initial T count " << running_t_count << std::endl;

    int lc = 0;
    int lb = 0;
    u64 time = 0;

    std::cout << "Loop count ";
    dT = join(dT, G, T, 0, &lb, &running_t_count, &time);

    lc++;
    while(true)
    {
      std::cout << "Loop count ";
      dT = join(dT, G, T, lc, &lb, &running_t_count, &time);

      if (lb == 1)  break;
      lc++;
    }

    delete dT;
    //delete T;

    char TCname[1024];
    sprintf(TCname, "%s_TC4", argv[1]);
    std::cout << "Filename " << TCname << std::endl;

    std::ofstream myfile;
    myfile.open (TCname);
    for (relation<2>::iter Tit(T, selectall); Tit.more(); Tit.advance())
      myfile << (*Tit)[0] << "\t" << (*Tit)[1] << "\n";
    myfile.close();

}



