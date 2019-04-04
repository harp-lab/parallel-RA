#include <chrono>
#include <math.h>
#include "btree.h"
#include "btree_relation.h"
#include "unordered_map"



uint32_t x=123456789, y=362436069, z=521288629;
inline uint32_t rnd()
{
    uint32_t t;
    x ^= x << 16;
    x ^= x >> 5;
    x ^= x << 1;

    t = x;
    x = y;
    y = z;
    z = t ^ x ^ y;

    return z;
}


void fatal(const char* const msg)
{
    std::cout << msg << std::endl;
    exit(0);
}


void testinsert(BTree<void>& t, u64 epochs, u64 max)
{
    for (uint64_t ep = 0; ep < epochs; ++ep)
        for (uint64_t j = 0; j < max; ++j)
        {
            const uint64_t i = max*ep + j;
            t.insert(i,(void*)1);

        }

    for (uint64_t ep = 0; ep < epochs; ++ep)
        for (uint64_t j = 0; j < max; ++j)
        {
            const uint64_t i = epochs*j + ep;
            t.insert(i,(void*)1);

        }

    for (uint64_t ep = 0; ep < epochs; ++ep)
        for (uint64_t j = 0; j < max; ++j)
        {
            const uint64_t i = epochs*j + ep;
            t.insert(i,(void*)1);

        }

    for (uint64_t ep = 0; ep < epochs; ++ep)
        for (uint64_t j = 0; j < max; ++j)
        {
            const uint64_t i = epochs*j + ep;
            t.insert(i,(void*)1);

        }
}

void testiterate(BTree<void>& t, u64 epochs, u64 max)
{
    u64 count = 0;
    for (BTree<void>::iter it(&t, -1); it.more(); it.advance())
        count++;

    if (count != epochs * max){
        std::cout << count << " " << epochs << " " << max << std::endl;
        fatal("BTree: Wrong element count");
    }
}

void testsearch(BTree<void>& t, u64 epochs, u64 max)
{
    u64 ccount = 0;
    u64 wcount = 0;

    for (uint64_t ep = 0; ep < epochs+2; ++ep)
        for (uint64_t j = 0; j < max; ++j)
        {
            const uint64_t i = max*ep + j;
            void* ans = t.search(i);
            if (ans == (void*)1)
                ccount++;
            else
                wcount++;
        }

    if (ccount != epochs * max)
        fatal("1 search error.");

    if (wcount != 2 * max)
        fatal("2 search error.");
}


void testinsert_map(std::unordered_map<u64, u64*>& t, u64 epochs, u64 max)
{
    for (uint64_t ep = 0; ep < epochs; ++ep)
        for (uint64_t j = 0; j < max; ++j)
        {
            const uint64_t i = max*ep + j;
            t.insert({i, (u64*)1});

        }

    for (uint64_t ep = 0; ep < epochs; ++ep)
        for (uint64_t j = 0; j < max; ++j)
        {
            const uint64_t i = epochs*j + ep;
            t.insert({i, (u64*)NULL});

        }

    for (uint64_t ep = 0; ep < epochs; ++ep)
        for (uint64_t j = 0; j < max; ++j)
        {
            const uint64_t i = epochs*j + ep;
            t.insert({i, (u64*)NULL});

        }

    for (uint64_t ep = 0; ep < epochs; ++ep)
        for (uint64_t j = 0; j < max; ++j)
        {
            const uint64_t i = epochs*j + ep;
            t.insert({i, (u64*)NULL});

        }
}

void testiterate_map(std::unordered_map<u64, u64*>& t, u64 epochs, u64 max)
{
    u64 count = 0;
    for ( auto it = t.begin(); it != t.end(); ++it )
        count++;

    if (count != epochs * max)
        fatal("Map: Wrong element count");
}

void testsearch_map(std::unordered_map<u64, u64*>& t, u64 epochs, u64 max)
{
    u64 ccount = 0;
    u64 wcount = 0;

    for (uint64_t ep = 0; ep < epochs+2; ++ep)
        for (uint64_t j = 0; j < max; ++j)
        {
            const uint64_t i = max*ep + j;
            auto ans = t.find(i);
            if (ans != t.end())
            {
                if (ans->second == (u64*)1)
                    ccount++;
            } else
                wcount++;
        }

    if (ccount != epochs * max)
    {
        std::cout <<  "Count: " << ccount << " " << epochs * max << std::endl;
        fatal("1 search error.");
    }

    if (wcount != 2 * max)
        fatal("2 search error.");
}


// Driver program to test above functions 
int main() 
{

    u64 epochs = 5;
    u64 max = 10000000;


    std::cout << "Unordered map" << std::endl;
    auto map_start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 1; i++)
    {
        std::unordered_map<u64, u64*> t;

        auto map_insert_start = std::chrono::high_resolution_clock::now();
        testinsert_map(t, epochs, max);
        auto map_insert_end = std::chrono::high_resolution_clock::now();


        auto map_iterate_start = std::chrono::high_resolution_clock::now();
        testiterate_map(t, epochs, max);
        auto map_iterate_end = std::chrono::high_resolution_clock::now();


        auto map_search_start = std::chrono::high_resolution_clock::now();
        testsearch_map(t, epochs, max);
        auto map_search_end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> map_insert_elapsed = map_insert_end - map_insert_start;
        std::chrono::duration<double> map_iterate_elapsed = map_iterate_end - map_iterate_start;
        std::chrono::duration<double> map_search_elapsed = map_search_end - map_search_start;

        std::cout << "Map insert " << map_insert_elapsed.count() << std::endl;
        std::cout << "Map iterate " << map_iterate_elapsed.count() << std::endl;
        std::cout << "Map search " << map_search_elapsed.count() << std::endl << std::endl;
    }
    auto map_end = std::chrono::high_resolution_clock::now();



    std::cout << "B tree" << std::endl;
    auto btree_start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 1; i++)
    {
        BTree<void> t(8);

        auto btree_insert_start = std::chrono::high_resolution_clock::now();
        testinsert(t, epochs, max);
        auto btree_insert_end = std::chrono::high_resolution_clock::now();


        auto btree_iterate_start = std::chrono::high_resolution_clock::now();
        testiterate(t, epochs, max);
        auto btree_iterate_end = std::chrono::high_resolution_clock::now();


        auto btree_search_start = std::chrono::high_resolution_clock::now();
        testsearch(t, epochs, max);
        auto btree_search_end = std::chrono::high_resolution_clock::now();


        std::chrono::duration<double> btree_insert_elapsed = btree_insert_end - btree_insert_start;
        std::chrono::duration<double> btree_iterate_elapsed = btree_iterate_end - btree_iterate_start;
        std::chrono::duration<double> btree_search_elapsed = btree_search_end - btree_search_start;

        std::cout << "Btree insert " << btree_insert_elapsed.count() << std::endl;
        std::cout << "Btree iterate " << btree_iterate_elapsed.count() << std::endl;
        std::cout << "Btree search " << btree_search_elapsed.count() << std::endl << std::endl;
    }
    auto btree_end = std::chrono::high_resolution_clock::now();


    std::cout << std::endl;
    std::cout << std::endl;

    std::chrono::duration<double> btree_elapsed = btree_end - btree_start;
    std::cout << "Btree " << btree_elapsed.count() << std::endl;

    std::cout << std::endl;

    std::chrono::duration<double> map_elapsed = map_end - map_start;
    std::cout << "Map " << map_elapsed.count() << std::endl;

    return 0;
}
