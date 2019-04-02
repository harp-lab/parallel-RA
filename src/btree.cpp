#include <chrono>
#include <math.h>
#include "btree.h"
#include "btree_relation.h"
#include "unordered_map"

/* The following program performs deletion on a B-Tree. It contains functions 
specific for deletion along with all the other functions provided in the 
previous articles on B-Trees. See https://www.geeksforgeeks.org/b-tree-set-1-introduction-2/ 
for previous article. 

The deletion function has been compartmentalized into 8 functions for ease 
of understanding and clarity 

The following functions are exclusive for deletion 
In class BTreeNode: 
        1) remove
        2) removeFromLeaf
        3) removeFromNonLeaf
        4) getPred
        5) getSucc
        6) borrowFromPrev
        7) borrowFromNext
        8) merge
        9) findKey

In class BTree: 
        1) remove

The removal of a key from a B-Tree is a fairly complicated process. The program handles 
all the 6 different cases that might arise while removing a key. 

Testing: The code has been tested using the B-Tree provided in the CLRS book( included 
in the main function ) along with other cases. 

Reference: CLRS3 - Chapter 18 - (499-502) 
It is advised to read the material in CLRS before taking a look at the code. */


uint64_t utime()
{
    return std::chrono::high_resolution_clock::now().time_since_epoch().count();
}


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
    //const uint64_t epochs = 5;
    //const uint64_t max = 100000;

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
            const uint64_t i = j*epochs + ep;
            void* ans = t.search(i);
            if (ans == (void*)1 && ep < epochs)
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
    //u64 v = 1;
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
    /*
    u64 count = 0;
    for ( auto it = t.begin(); it != t.end(); ++it )
    {
       if (it->first < 0 || it->first > ecount)
           fatal("2 Bad value during iteration.");
       count++;
    }

    if (count != ecount)
        fatal("Wrong element count");
    */

    u64 ccount = 0;
    u64 wcount = 0;

    for (uint64_t ep = 0; ep < epochs+2; ++ep)
        for (uint64_t j = 0; j < max; ++j)
        {
            const uint64_t i = epochs*j + ep;
            u64* ans = t[i];
            if (ans == NULL && ep < epochs)
                ccount++;
            else
                wcount++;
        }

    if (ccount != epochs * max)
        fatal("1 search error.");

    if (wcount != 2 * max)
        fatal("2 search error.");
}


void testrandom(BTree<void>& t)
{
    for (uint64_t i = 0; i < 1000000; ++i)
    {
        const uint64_t n = (((uint64_t)rnd()) << 32) | rnd();
        t.insert(n,(void*)1);
    }

    for (uint64_t i = 100000; i < 500000; ++i)
    {
        const uint64_t n = (((uint64_t)rnd()) << 32) | rnd();
        //const uint64_t n = ((uint64_t)rnd()) % 100;
        t.remove(n);
    }

    for (uint64_t i = 0; i < 700000; ++i)
    {
        const uint64_t n = (((uint64_t)rnd()) << 32) | rnd();
        t.search(n);
    }
}

void test(BTree<void>& t)
{
    FILE *fp_in;
    fp_in = fopen("111", "r");
    int element1 = 10, element2 = 10;
    int rc = 1142;

    for (int i = 0; i < rc; i++)
    {
        if (fscanf(fp_in,"%d\t%d\n", &element1, &element2) !=2 )
            break;
        //fscanf(fp_in,"%d\t%d\n", &element1, &element2);

        bool x = t.insert((u64)element2,(void*)1);

        std::cout << element2 << "\t" << x << std::endl;
    }
    fclose(fp_in);

    for (BTree<void>::iter it(&t, -1); it.more(); it.advance())
      std::cout << it.getkey() << std::endl;
}


void testseq(BTree<void>& t)
{
    const uint64_t epochs = 5;
    const uint64_t max = 100000;

    for (uint64_t ep = 0; ep < epochs; ++ep)
        for (uint64_t j = 0; j < max; ++j)
        {
            const uint64_t i = ep*0x300000 + j;
            t.insert(i,(void*)1);

        }



    for (BTree<void>::iter it(&t, -1); it.more(); it.advance())
    {

        if (it.getval() != (void*)1)
            fatal("Bad value during iteration.");
        for (uint64_t ep = 0; ep < epochs; ++ep)
            if (it.getkey() >= ep*0x300000 && it.getkey() < ep*0x300000+max)
                goto continuelab;
        fatal("Key not in range during iteration");
continuelab:;
    }


    /*

    for (uint64_t ep = 1; ep < epochs; ++ep)
        for (uint64_t j = 0; j < max; ++j)
        {
            const uint64_t i = ep*0x300000 + j;
            t.remove(i);
        }
    */

    for (uint64_t ep = 0; ep < epochs+2; ++ep)
        for (uint64_t j = 0; j < max; ++j)
        {
            const uint64_t i = ep*0x300000 + j;
            void* ans = t.search(i);
            if (ans == 0 && ep == 0)
                fatal("Should have found key.");
            else if (ans != 0 && ep > 0)
                fatal("Should not have found key.");
        }


    //printf("count = %d\n", (int)t.count());

    if (t.count() != max)
        fatal("Wrong key-count.");


}

// Driver program to test above functions 
int main() 
{

    u64 epochs = 5;
    u64 max = 10000000;

    std::cout << "B tree" << std::endl;
    for (int i = 1; i < 5; i++)
    {
        BTree<void> t(8);

        uint64_t start = utime();
        testinsert(t, epochs, max);
        uint64_t end = utime();

        std::cout << "Insert: " << max * epochs << " in " << ((end - start)/1000000) << std::endl;

        testiterate(t, epochs, max);
        uint64_t end2 = utime();
        std::cout << "Iterate: " << ((end2 - end)/1000000) << std::endl;

        testsearch(t, epochs, max);
        uint64_t end3 = utime();
        std::cout << "Iterate and search: " << ((end3 - end2)/1000000) << std::endl;
        std::cout << "Total time: " << (end3 - start)/1000000 << std::endl << std::endl;
    }

    std::cout << "Unordered map" << std::endl;
    for (int i = 1; i < 5; i++)
    {
        std::unordered_map<u64, u64*> t;

        uint64_t start = utime();
        testinsert_map(t, epochs, max);
        uint64_t end = utime();

        std::cout << "Insert: " << max * epochs << " in " << ((end - start)/1000000) << std::endl;

        testiterate_map(t, epochs, max);
        uint64_t end2 = utime();
        std::cout << "Iterate: " << ((end2 - end)/1000000) << std::endl;

        testsearch_map(t, epochs, max);
        uint64_t end3 = utime();
        std::cout << "Iterate and search: " << ((end3 - end2)/1000000) << std::endl;

        std::cout << "Total time: " << (end3 - start)/1000000 << std::endl << std::endl;
    }


    /*
    u64 start = utime();
    for (uint64_t i = 0; i < 1; ++i)
    {
        BTree<void> t(8);
        //testseq(t);
        test(t);
    }
    std::cout << ((utime() - start)/1000000) << std::endl;
    */

    return 0;
}
