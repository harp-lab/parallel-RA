#include <chrono>
#include "btree.h"
#include "btree_relation.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include <unordered_map>
#include <unordered_set>



void fatal(const char* const msg)
{
    std::cout << msg << std::endl;
    exit(0);
}


void fileIO(const char* filename, tuple<2> * tup, int rc)
{
    FILE *fp_in  = fopen(filename, "r");

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

}



void testUmap(tuple<2> * tup, u64 rc)
{
    //auto insert_start = std::chrono::high_resolution_clock::now();
    std::unordered_map<u64, std::unordered_set<u64>> G;
    u64 row_count = 0;
    for (u64 i = 0; i < rc; i++)
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

    //if (row_count != rc)
    //    std::cout << "Error !!!!" << row_count << " " << rc << std::endl;

    //auto insert_finish = std::chrono::high_resolution_clock::now();


    /*
    row_count = 0;
    std::unordered_set<u64> k;
    for ( auto local_it = G.begin(); local_it!= G.end(); ++local_it )
    {
        k = local_it->second;
        for (auto it2 = k.begin(); it2 != k.end(); it2++)
            row_count++;
    }

    if (row_count != rc)
        std::cout << "Error !!!!" << row_count << " " << rc << std::endl;
    */

    //auto iterate_finish = std::chrono::high_resolution_clock::now();

    //std::chrono::duration<double> insert_elapsed = insert_finish - insert_start;
    //std::chrono::duration<double> iterate_elapsed = iterate_finish - insert_finish;

    //std::cout << "Map insert " << insert_elapsed.count() << std::endl;
    //std::cout << "Map iterate " << iterate_elapsed.count() << std::endl;
}




void testBtree(tuple<2> * tup, u64 rc)
{
    //auto insert_start = std::chrono::high_resolution_clock::now();
    u64 row_count = 0;
    relation<2> rel;

    for (u64 i = 0; i < rc; i++)
    {
        if (rel.insert(tup[i]))
          row_count++;
    }

    //if (row_count != rc)
    //    std::cout << "Error !!!!" << row_count << " " << rc << std::endl;

    //auto insert_finish = std::chrono::high_resolution_clock::now();

    /*
    tuple<2> t1;
    t1[0] = -1; t1[1] = -1;
    tuple<2> selectall(t1);
    row_count = 0;
    for (relation<2>::iter it(rel, selectall); it.more(); it.advance())
        row_count++;

    if (row_count != rc)
        std::cout << "Error !!!!" << row_count << " " << rc << std::endl;
    */

    //auto iterate_finish = std::chrono::high_resolution_clock::now();

    //std::chrono::duration<double> insert_elapsed = insert_finish - insert_start;
    //std::chrono::duration<double> iterate_elapsed = iterate_finish - insert_finish;

    //std::cout << "Btree insert " << insert_elapsed.count() << std::endl;
    //std::cout << "Btree iterate " << iterate_elapsed.count() << std::endl;
}



// Driver program to test above functions
int main()
{
    tuple<2> * tup = new tuple<2>[16411677];
    fileIO("24_T", tup, 16411677);

    auto map_start = std::chrono::high_resolution_clock::now();
    testUmap(tup, 16411677);
    auto map_finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> map_elapsed = map_finish - map_start;
    std::cout << "16411677 Map insert " << map_elapsed.count() << std::endl;

    auto btree_start = std::chrono::high_resolution_clock::now();
    testBtree(tup, 16411677);
    auto btree_finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> btree_elapsed = btree_finish - btree_start;
    std::cout << "16411677 Btree insert " << btree_elapsed.count() << std::endl;




    std::cout << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << std::endl;


    delete[] tup;


    return 0;
}
