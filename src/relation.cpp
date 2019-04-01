#include <chrono>
#include "btree.h"
#include "btree_relation.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>

uint64_t utime()
{
    return std::chrono::high_resolution_clock::now().time_since_epoch().count();
}



void fatal(const char* const msg)
{
    std::cout << msg << std::endl;
    exit(0);
}




void testseq()
{
    relation<2> rel;
    tuple<2> t1;

    FILE *fp_in;
    //fp_in = fopen("link.facts_412148", "r");
    fp_in = fopen("X1", "r");
    assert (fp_in != NULL);
    int row_count = 0;
    int element1 = 10, element2 = 10;
    int rc = 412148;

    for (int i = 0; i < rc; i++)
    {
        if (fscanf(fp_in,"%d\t%d\n", &element1, &element2) !=2 )
            break;
        //fscanf(fp_in,"%d\t%d\n", &element1, &element2);

        t1[0] = (u64)element1;
        t1[1] = (u64)element2;

        //std::cout << "[EXE] " << element1 << ", " << element2 << std::endl;
        //bool x = rel.insert(t1);
        rel.insert(t1);

        //if (x == true && t1[0] == 11 && t1[1] == 134)
        //    std::cout << "Shout: " << t1[0] << "\t" << t1[1] << std::endl;

        row_count++;
    }
    fclose(fp_in);


    relation<2> * dT = new relation<2>;
    t1[0] = -1; t1[1] = -1;
    tuple<2> selectall(t1);
    for (relation<2>::iter it(rel, selectall); it.more(); it.advance())
    {
        tuple<2> t1;
        t1[0] = (*it)[0];
        t1[1] = (*it)[1];
        if (dT->insert(t1) == true)
            ;
        else
            std::cout << "Shout" << std::endl;
      //std::cout << (*it)[0] << "\t" << (*it)[1] << std::endl;
    }
}

// Driver program to test above functions
int main()
{
    u64 start = utime();
    testseq();
    std::cout << ((utime() - start)/1000000) << std::endl;

    return 0;
}
