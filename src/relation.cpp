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
#include "relation.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "btree/btree_map.h"

#include <unordered_map>
#include <unordered_set>


void testUmap(u64 * tup, u64 rc)
{
    std::unordered_map<u64, std::unordered_set<u64>*> *G = new std::unordered_map<u64, std::unordered_set<u64>*>;
    u64 row_count = 0;
    for (u64 i = 0; i < rc * 2; i=i+2)
    {
        auto it = G->find(tup[i]);
        if( it != G->end() ) {
            auto it2 = (it->second)->find(tup[i + 1]);
            if( it2 != (it->second)->end() ) {
                ;
            }
            else{
                (it->second)->insert(tup[i + 1]);
                (*G)[tup[i]] = it->second;
                row_count++;
            }
        }
        else {
            std::unordered_set<u64>* k = new std::unordered_set<u64>;
            k->insert(tup[i + 1]);
            G->insert(std::make_pair(tup[i],k));
            row_count++;
        }
    }

    if (row_count != rc)
        std::cout << "Error !!!!" << row_count << " " << rc << std::endl;

    std::unordered_map<u64, std::unordered_set<u64>*>::iterator i1 = G->begin();
    for(; i1 != G->end(); i1++)
        delete (i1->second);

    delete G;
}


void testBtree(u64 * tup, u64 rc)
{
    u64 row_count = 0;
    relation<2> *rel = new relation<2>;
    tuple<2> t1;
    for (u64 i = 0; i < rc * 2; i=i+2)
    {
        t1[0] = tup[i]; t1[1] = tup[i+1];
        if (rel->insert(t1))
          row_count++;
    }

    if (row_count != rc)
        std::cout << "Error !!!!" << row_count << " " << rc << std::endl;

    delete rel;
}



void testGoogleBtree(u64 * tup, u64 rc)
{
    typedef btree::btree_map<u64, u64> Relation0Map;
    typedef btree::btree_map<u64, btree::btree_map<u64, u64>* > Relation1Map;


    Relation1Map *map1 = new Relation1Map;

    u64 row_count = 0;
    for (u64 i = 0; i < rc * 2; i=i+2)
    {
        Relation1Map::iterator it = map1->find(tup[i]);
        if ( it == map1->end() )
        {
            Relation0Map *map0 = new Relation0Map;
            map0->insert(std::make_pair(tup[i + 1], 0));
            map1->insert(std::make_pair(tup[i], map0));
            row_count++;
        }
        else
        {
          auto it2 = (it->second)->find(tup[i + 1]);
          if( it2 != (it->second)->end() ) {
              ;
          }
          else{

              Relation0Map *map0 = (it->second);
              map0->insert(std::make_pair(tup[i + 1], 0));
              map1->insert(std::make_pair(tup[i], map0));
              row_count++;
          }
        }
    }

    if (row_count != rc)
        std::cout << "Error !!!!" << row_count << " " << rc << std::endl;


    Relation1Map::iterator i1 = map1->begin();
    for(; i1 != map1->end(); i1++)
        delete (i1->second);

    delete map1;
}


void check(u32 index, u64 rc, u64* read_buffer)
{

    //  std::cout << std::endl;
    //std::cout << index << " [O] ------------------------------------------- " << rc << " -------------------------------------------" << std::endl;


    /*
    auto btree_start = std::chrono::high_resolution_clock::now();
    for (u32 i = 0; i < iteration_count; i++)
    {
        //auto btree_start_i = std::chrono::high_resolution_clock::now();
        testBtree(read_buffer, rc);
        //auto btree_finish_i = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> btree_elapsed_i = btree_finish_i - btree_start_i;
        //std::cout << "[" << rc << "] [" << i <<"] Btree insert " << btree_elapsed_i.count() << std::endl;
    }
    auto btree_finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> btree_elapsed = btree_finish - btree_start;
    std::cout << "[" << rc << "] Btree " << btree_elapsed.count() / iteration_count << std::endl;
    */


    u32 iteration_count = 5;
    auto google_start = std::chrono::high_resolution_clock::now();
    for (u32 i = 0; i < iteration_count; i++)
    {
        //auto google_start_i = std::chrono::high_resolution_clock::now();
        testGoogleBtree(read_buffer, rc);
        //auto google_finish_i = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> google_elapsed_i = google_finish_i - google_start_i;
        //std::cout << "[" << rc << "] [" << i <<"] Google insert " << google_elapsed_i.count() << std::endl;
    }
    auto google_finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> google_elapsed = google_finish - google_start;
    std::cout << rc << "\t" << "Google" << "\t" << google_elapsed.count() / iteration_count << std::endl;



    auto map_start = std::chrono::high_resolution_clock::now();
    for (u32 i = 0; i < iteration_count; i++)
    {
        //auto map_start_i = std::chrono::high_resolution_clock::now();
        testUmap(read_buffer, rc);
        //auto map_finish_i = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> map_elapsed_i = map_finish_i - map_start_i;
        //std::cout << "[" << rc << "] [" << i <<"] Map insert " << map_elapsed_i.count() << std::endl;
    }
    auto map_finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> map_elapsed = map_finish - map_start;
    //std::cout << "[" << rc << "] Map " << map_elapsed.count() / iteration_count << std::endl;
    std::cout << rc << "\t" << "Table" << "\t" << map_elapsed.count() / iteration_count << std::endl;



    std::cout << std::endl;
    std::cout << std::endl;


}


// Driver program to test above functions
int main()
{

    u64 rowC = 40451631;
    int fp = open("data.raw", O_RDONLY);
    u64 * read_buffer = new u64[rowC * 2];
    u32 rb_size = pread(fp, read_buffer, rowC * 2 * sizeof(u64), 0);
    if (rb_size != rowC * 2 * sizeof(u64))
      std::cout << "Error !!!!" << rb_size << " " << rowC * 2 * sizeof(u64) << std::endl;
    close(fp);

    //check(0, rowC, read_buffer);

    u32 c = 0;
    for (u64 rc = 0; rc <= rowC; rc = rc +  400000)
    {
        check(c, rc, read_buffer);
        c++;
    }
    check(c,136024430 , read_buffer);

    delete[] read_buffer;


    return 0;
}
