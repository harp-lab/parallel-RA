#include <sys/stat.h>
#include <errno.h>
#include <limits.h>
#include <arpa/inet.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
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


/*
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
*/



void testUmap(u64 * tup, u64 rc)
{
    //auto insert_start = std::chrono::high_resolution_clock::now();
    std::unordered_map<u64, std::unordered_set<u64>> G;
    u64 row_count = 0;
    for (u64 i = 0; i < rc * 2; i=i+2)
    {
        auto it = G.find(tup[i]);
        if( it != G.end() ) {
            auto it2 = (it->second).find(tup[i + 1]);
            if( it2 != (it->second).end() ) {
                ;
            }
            else{
                (it->second).insert(tup[i + 1]);
                G[tup[i]] = it->second;
                row_count++;
            }
        }
        else {
            std::unordered_set<u64> k;
            k.insert(tup[i + 1]);
            G.insert(std::make_pair(tup[i],k));
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




void testBtree(u64 * tup, u64 rc)
{
    //auto insert_start = std::chrono::high_resolution_clock::now();
    u64 row_count = 0;
    relation<2> rel;
    tuple<2> t1;
    for (u64 i = 0; i < rc * 2; i=i+2)
    {
        t1[0] = tup[i]; t1[1] = tup[i+1];
        if (rel.insert(t1))
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


void check(u64 rc, u64* read_buffer)
{

    std::cout << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << std::endl;
    std::cout << "Tuple count: " << rc << std::endl;
    std::cout << std::endl;

    auto map_start = std::chrono::high_resolution_clock::now();
    testUmap(read_buffer, rc);
    auto map_finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> map_elapsed = map_finish - map_start;
    std::cout << "Map insert " << map_elapsed.count() << std::endl;

    auto btree_start = std::chrono::high_resolution_clock::now();
    testBtree(read_buffer, rc);
    auto btree_finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> btree_elapsed = btree_finish - btree_start;
    std::cout << "Btree insert " << btree_elapsed.count() << std::endl;


    map_start = std::chrono::high_resolution_clock::now();
    testUmap(read_buffer, rc);
    map_finish = std::chrono::high_resolution_clock::now();
    map_elapsed = map_finish - map_start;
    std::cout << "Map insert " << map_elapsed.count() << std::endl;

    btree_start = std::chrono::high_resolution_clock::now();
    testBtree(read_buffer, rc);
    btree_finish = std::chrono::high_resolution_clock::now();
    btree_elapsed = btree_finish - btree_start;
    std::cout << "Btree insert " << btree_elapsed.count() << std::endl;


    map_start = std::chrono::high_resolution_clock::now();
    testUmap(read_buffer, rc);
    map_finish = std::chrono::high_resolution_clock::now();
    map_elapsed = map_finish - map_start;
    std::cout << "Map insert " << map_elapsed.count() << std::endl;

    btree_start = std::chrono::high_resolution_clock::now();
    testBtree(read_buffer, rc);
    btree_finish = std::chrono::high_resolution_clock::now();
    btree_elapsed = btree_finish - btree_start;
    std::cout << "Btree insert " << btree_elapsed.count() << std::endl;



    map_start = std::chrono::high_resolution_clock::now();
    testUmap(read_buffer, rc);
    map_finish = std::chrono::high_resolution_clock::now();
    map_elapsed = map_finish - map_start;
    std::cout << "Map insert " << map_elapsed.count() << std::endl;

    btree_start = std::chrono::high_resolution_clock::now();
    testBtree(read_buffer, rc);
    btree_finish = std::chrono::high_resolution_clock::now();
    btree_elapsed = btree_finish - btree_start;
    std::cout << "Btree insert " << btree_elapsed.count() << std::endl;


    map_start = std::chrono::high_resolution_clock::now();
    testUmap(read_buffer, rc);
    map_finish = std::chrono::high_resolution_clock::now();
    map_elapsed = map_finish - map_start;
    std::cout << "Map insert " << map_elapsed.count() << std::endl;

    btree_start = std::chrono::high_resolution_clock::now();
    testBtree(read_buffer, rc);
    btree_finish = std::chrono::high_resolution_clock::now();
    btree_elapsed = btree_finish - btree_start;
    std::cout << "Btree insert " << btree_elapsed.count() << std::endl;


    map_start = std::chrono::high_resolution_clock::now();
    testUmap(read_buffer, rc);
    map_finish = std::chrono::high_resolution_clock::now();
    map_elapsed = map_finish - map_start;
    std::cout << "Map insert " << map_elapsed.count() << std::endl;

    btree_start = std::chrono::high_resolution_clock::now();
    testBtree(read_buffer, rc);
    btree_finish = std::chrono::high_resolution_clock::now();
    btree_elapsed = btree_finish - btree_start;
    std::cout << "Btree insert " << btree_elapsed.count() << std::endl;



    std::cout << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << std::endl;

}


// Driver program to test above functions
int main()
{

    u64 rc = 32823354/2;
    int fp = open("24_T", O_RDONLY);
    u64 * read_buffer = new u64[rc * 2];
    u32 rb_size = pread(fp, read_buffer, rc * 2 * sizeof(u64), 0);
    if (rb_size != rc * 2 * sizeof(u64))
    {
        std::cout << "Error !!!!" << std::endl;
    }
    close(fp);


    rc = 1000;
    check (rc, read_buffer);

    rc = 5000;
    check (rc, read_buffer);

    rc = 10000;
    check (rc, read_buffer);

    rc = 50000;
    check (rc, read_buffer);

    rc = 100000;
    check (rc, read_buffer);

    rc = 500000;
    check (rc, read_buffer);

    rc = 1000000;
    check (rc, read_buffer);


    delete[] read_buffer;


    return 0;
}
