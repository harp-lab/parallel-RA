#include <chrono>
#include "btree.h"

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


void testrandom(BTree& t)
{
    for (uint64_t i = 0; i < 500000; ++i)
    {
        const uint64_t n = (((uint64_t)rnd()) << 32) | rnd();
        t.insert(n,(BTreeNode*)1);
    }
    for (uint64_t i = 100000; i < 500000; ++i)
    {
        const uint64_t n = (((uint64_t)rnd()) << 32) | rnd();
        t.remove(n);
    }
    for (uint64_t i = 0; i < 700000; ++i)
    {
        const uint64_t n = (((uint64_t)rnd()) << 32) | rnd();
        t.search(n);
    }
}


void testseq(BTree& t)
{
    const uint64_t epochs = 5;
    const uint64_t max = 100000;
    for (uint64_t ep = 0; ep < epochs; ++ep)
        for (uint64_t j = 0; j < max; ++j)
        {
            const uint64_t i = ep*0x300000 + j;
            t.insert(i,(BTreeNode*)1);

        }


    for (BTree::iter it(&t); it.more(); it.advance())
    {

        if (it.getval() != (void*)1)
            fatal("Bad value during iteration.");
        for (uint64_t ep = 0; ep < epochs; ++ep)
            if (it.getkey() >= ep*0x300000 && it.getkey() < ep*0x300000+max)
                goto continuelab;
        fatal("Key not in range during iteration");
continuelab:;
    }


    for (uint64_t ep = 1; ep < epochs; ++ep)
        for (uint64_t j = 0; j < max; ++j)
        {
            const uint64_t i = ep*0x300000 + j;
            t.remove(i);
        }

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
    uint64_t start = utime();
    for (uint64_t i = 0; i < 3; ++i)
    {
        BTree t(8);
        testrandom(t);
    }
    std::cout << ((utime() - start)/1000000) << std::endl;

    start = utime();
    for (uint64_t i = 0; i < 1; ++i)
    {
        BTree t(8);
        testseq(t);
    }
    std::cout << ((utime() - start)/1000000) << std::endl;
    return 0;
}
