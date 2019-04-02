#include <chrono>
#include "btree.h"
#include "btree_relation.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include <unordered_map>
#include <unordered_set>

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
    tuple<2> * tup = new tuple<2>[40451631];

    FILE *fp_in  = fopen("d_40451631", "r");

    int row_count = 0;
    int rc = 40451631;
    int element1 = 0;
    int element2 = 0;

    for (int i = 0; i < rc; i++)
    {
        if (fscanf(fp_in,"%d\t%d\n", &element1, &element2) !=2 )
            break;

        tup[i][0] = (u64)element1;
        tup[i][1] = (u64)element2;
    }
    fclose(fp_in);


    u64 start = utime();
    for (int i = 0; i < rc; i++)
    {
        if (rel.insert(tup[i]))
          row_count++;
    }
    std::cout << "Btree insert: " << ((utime() - start)/1000000) << std::endl;

    if (row_count != rc)
        std::cout << "Error !!!!" << row_count << " " << rc << std::endl;

    tuple<2> t1;
    t1[0] = -1; t1[1] = -1;
    tuple<2> selectall(t1);
    row_count = 0;
    start = utime();
    for (relation<2>::iter it(rel, selectall); it.more(); it.advance())
        row_count++;
    std::cout << "Btree iterator: " << ((utime() - start)/1000000) << std::endl;

    if (row_count != rc)
        std::cout << "Error !!!!" << row_count << " " << rc << std::endl;


    std::unordered_map<u64, std::unordered_set<u64>> G;
    row_count = 0;
    start = utime();
    for (int i = 0; i < rc; i++)
    {
        auto it = G.find(tup[i][0]);
        if( it != G.end() ) {
            auto it2 = (it->second).find(tup[i][1]);
            if( it2 != (it->second).end() ) {
                ;
            }
            else{
                (it->second).insert(tup[i][1]);
                G[tup[i][0]] = it->second;
                row_count++;
            }
        }
        else {
            std::unordered_set<u64> k;
            k.insert(tup[i][1]);
            G.insert(std::make_pair(tup[i][0],k));
            row_count++;
        }
    }

    if (row_count != rc)
        std::cout << "Error !!!!" << row_count << " " << rc << std::endl;

    std::cout << "Map insert " << ((utime() - start)/1000000) << std::endl;

    row_count = 0;
    start = utime();
    std::unordered_set<u64> k;
    for ( auto local_it = G.begin(); local_it!= G.end(); ++local_it )
    {
        k = local_it->second;
        for (auto it2 = k.begin(); it2 != k.end(); it2++)
            row_count++;
    }
    std::cout << "Map iterator " << ((utime() - start)/1000000) << std::endl;

    if (row_count != rc)
        std::cout << "Error !!!!" << row_count << " " << rc << std::endl;

    delete[] tup;

#if 0
    tuple<2> t1;
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
#endif

}

// Driver program to test above functions
int main()
{

    testseq();


    return 0;
}
