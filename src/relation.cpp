#include <chrono>
#include "btree.h"
#include "relation.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>




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




void testseq()
{
    relation<2> rel;
    tuple<2> t1;

    FILE *fp_in;
    fp_in = fopen("TCS", "r");
    int row_count = 100;
    long long int element1 = 0, element2 = 0;
    //int rc = 1711;
    int rc = 46;

    for (int i = 0; i < rc; i++)
    {
        assert(fscanf(fp_in,"%lld\t%lld\n", &element1, &element2) == 2);

        t1[0] = (u64)element1;
        t1[1] = (u64)element2;

        //std::cout << "[EXE] " << t1[0] << ", " << t1[1] << std::endl;
        rel.insert(t1);


        row_count++;
    }
    fclose(fp_in);

    /*
    t1[0] = 2;
    t1[1] = 1;
    rel.insert(t1);

    t1[0] = 2;
    t1[1] = 2;
    rel.insert(t1);

    t1[0] = 2;
    t1[1] = 3;
    rel.insert(t1);

    t1[0] = 3;
    t1[1] = 4;
    rel.insert(t1);

    t1[0] = 3;
    t1[1] = 6;
    rel.insert(t1);

    t1[0] = 3;
    t1[1] = 8;
    rel.insert(t1);

    t1[0] = 4;
    t1[1] = 19;
    rel.insert(t1);

    t1[0] = 4;
    t1[1] = 2;
    rel.insert(t1);

    t1[0] = 4;
    t1[1] = 20;
    rel.insert(t1);

    t1[0] = 5;
    t1[1] = 19;
    rel.insert(t1);

    t1[0] = 5;
    t1[1] = 21;
    rel.insert(t1);

    t1[0] = 5;
    t1[1] = 9;
    rel.insert(t1);

    t1[0] = 6;
    t1[1] = 17;
    rel.insert(t1);

    t1[0] = 6;
    t1[1] = 21;
    rel.insert(t1);

    t1[0] = 6;
    t1[1] = 22;
    rel.insert(t1);

    t1[0] = 7;
    t1[1] = 24;
    rel.insert(t1);

    t1[0] = 7;
    t1[1] = 25;
    rel.insert(t1);

    t1[0] = 7;
    t1[1] = 5;
    rel.insert(t1);

    t1[0] = 8;
    t1[1] = 10;
    rel.insert(t1);

    t1[0] = 8;
    t1[1] = 10;
    rel.insert(t1);

    t1[0] = 8;
    t1[1] = 26;
    rel.insert(t1);

    t1[0] = 8;
    t1[1] = 26;
    rel.insert(t1);

    t1[0] = 8;
    t1[1] = 27;
    rel.insert(t1);

    t1[0] = 8;
    t1[1] = 27;
    rel.insert(t1);

    t1[0] = 9;
    t1[1] = 11;
    rel.insert(t1);

    t1[0] = 9;
    t1[1] = 12;
    rel.insert(t1);

    t1[0] = 9;
    t1[1] = 4;
    rel.insert(t1);

    t1[0] = 10;
    t1[1] = 1;
    rel.insert(t1);

    t1[0] = 10;
    t1[1] = 17;
    rel.insert(t1);

    t1[0] = 10;
    t1[1] = 21;
    rel.insert(t1);

    t1[0] = 10;
    t1[1] = 22;
    rel.insert(t1);

    t1[0] = 10;
    t1[1] = 6;
    rel.insert(t1);

    t1[0] = 10;
    t1[1] = 7;
    rel.insert(t1);

    t1[0] = 11;
    t1[1] = 1;
    rel.insert(t1);

    t1[0] = 11;
    t1[1] = 1;
    rel.insert(t1);

    t1[0] = 11;
    t1[1] = 2;
    rel.insert(t1);

    t1[0] = 11;
    t1[1] = 20;
    rel.insert(t1);

    t1[0] = 11;
    t1[1] = 24;
    rel.insert(t1);

    t1[0] = 11;
    t1[1] = 3;
    rel.insert(t1);

    t1[0] = 11;
    t1[1] = 30;
    rel.insert(t1);

    t1[0] = 11;
    t1[1] = 6;
    rel.insert(t1);

    t1[0] = 11;
    t1[1] = 7;
    rel.insert(t1);

    t1[0] = 12;
    t1[1] = 12;
    rel.insert(t1);

    t1[0] = 12;
    t1[1] = 14;
    rel.insert(t1);

    t1[0] = 12;
    t1[1] = 15;
    rel.insert(t1);

    t1[0] = 12;
    t1[1] = 15;
    rel.insert(t1);

    t1[0] = 12;
    t1[1] = 18;
    rel.insert(t1);

    t1[0] = 12;
    t1[1] = 4;
    rel.insert(t1);

    t1[0] = 12;
    t1[1] = 5;
    rel.insert(t1);

    t1[0] = 12;
    t1[1] = 6;
    rel.insert(t1);

    t1[0] = 12;
    t1[1] = 8;
    rel.insert(t1);

    t1[0] = 13;
    t1[1] = 11;
    rel.insert(t1);

    t1[0] = 13;
    t1[1] = 12;
    rel.insert(t1);

    t1[0] = 13;
    t1[1] = 13;
    rel.insert(t1);

    t1[0] = 13;
    t1[1] = 14;
    rel.insert(t1);

    t1[0] = 13;
    t1[1] = 19;
    rel.insert(t1);

    t1[0] = 13;
    t1[1] = 21;
    rel.insert(t1);

    t1[0] = 13;
    t1[1] = 3;
    rel.insert(t1);

    t1[0] = 13;
    t1[1] = 4;
    rel.insert(t1);

    t1[0] = 13;
    t1[1] = 9;
    rel.insert(t1);

    t1[0] = 14;
    t1[1] = 10;
    rel.insert(t1);

    t1[0] = 14;
    t1[1] = 13;
    rel.insert(t1);

    t1[0] = 14;
    t1[1] = 14;
    rel.insert(t1);

    t1[0] = 14;
    t1[1] = 16;
    rel.insert(t1);

    t1[0] = 14;
    t1[1] = 17;
    rel.insert(t1);

    t1[0] = 14;
    t1[1] = 25;
    rel.insert(t1);

    t1[0] = 14;
    t1[1] = 26;
    rel.insert(t1);

    t1[0] = 14;
    t1[1] = 3;
    rel.insert(t1);

    t1[0] = 14;
    t1[1] = 8;
    rel.insert(t1);

    t1[0] = 15;
    t1[1] = 10;
    rel.insert(t1);

    t1[0] = 15;
    t1[1] = 13;
    rel.insert(t1);

    t1[0] = 15;
    t1[1] = 14;
    rel.insert(t1);

    t1[0] = 15;
    t1[1] = 15;
    rel.insert(t1);

    t1[0] = 15;
    t1[1] = 16;
    rel.insert(t1);

    t1[0] = 15;
    t1[1] = 26;
    rel.insert(t1);

    t1[0] = 15;
    t1[1] = 27;
    rel.insert(t1);

    t1[0] = 15;
    t1[1] = 5;
    rel.insert(t1);

    t1[0] = 15;
    t1[1] = 7;
    rel.insert(t1);

    t1[0] = 16;
    t1[1] = 10;
    rel.insert(t1);

    t1[0] = 16;
    t1[1] = 13;
    rel.insert(t1);

    t1[0] = 16;
    t1[1] = 16;
    rel.insert(t1);

    t1[0] = 16;
    t1[1] = 16;
    rel.insert(t1);

    t1[0] = 16;
    t1[1] = 17;
    rel.insert(t1);

    t1[0] = 16;
    t1[1] = 18;
    rel.insert(t1);

    t1[0] = 16;
    t1[1] = 28;
    rel.insert(t1);

    t1[0] = 16;
    t1[1] = 29;
    rel.insert(t1);

    t1[0] = 16;
    t1[1] = 7;
    rel.insert(t1);

    t1[0] = 17;
    t1[1] = 11;
    rel.insert(t1);

    t1[0] = 17;
    t1[1] = 12;
    rel.insert(t1);

    t1[0] = 17;
    t1[1] = 15;
    rel.insert(t1);

    t1[0] = 17;
    t1[1] = 18;
    rel.insert(t1);

    t1[0] = 17;
    t1[1] = 19;
    rel.insert(t1);

    t1[0] = 17;
    t1[1] = 2;
    rel.insert(t1);

    t1[0] = 17;
    t1[1] = 20;
    rel.insert(t1);

    t1[0] = 17;
    t1[1] = 23;
    rel.insert(t1);

    t1[0] = 17;
    t1[1] = 28;
    rel.insert(t1);


    t1[0] = 18;
    t1[1] = 11;
    rel.insert(t1);

    t1[0] = 18;
    t1[1] = 23;
    rel.insert(t1);

    t1[0] = 18;
    t1[1] = 28;
    rel.insert(t1);

    t1[0] = 18;
    t1[1] = 9;
    rel.insert(t1);
    */

    //std::cout << "inserted rows: " << row_count << std::endl;

    t1[0] = -1; t1[1] = -1;
    tuple<2> selectall(t1);
    /*
    relation<2>::iter it(rel, selectall);
    std::cout << "[A] " << (*it)[0] << "\t" << (*it)[1] << std::endl;
    it.advance();
    std::cout << "[B] " << (*it)[0] << "\t" << (*it)[1] << std::endl;
    it.advance();
    std::cout << "[C] " << (*it)[0] << "\t" << (*it)[1] << std::endl;
    it.advance();
    std::cout << "[D] " << (*it)[0] << "\t" << (*it)[1] << std::endl;
    it.advance();
    std::cout << "[E] " << (*it)[0] << "\t" << (*it)[1] << std::endl;
    it.advance();
    std::cout << "[F] " << (*it)[0] << "\t" << (*it)[1] << std::endl;
    */

    for (relation<2>::iter it(rel, selectall); it.more(); it.advance())
      std::cout << (*it)[0] << "\t" << (*it)[1] << std::endl;

    /*
    const uint64_t epochs = 5;
    const uint64_t max = 10000;
    for (uint64_t ep = 0; ep < epochs; ++ep)
        for (uint64_t j = 0; j < max; ++j)
        {
            const uint64_t i = ep*0x300000 + j;

            t1[0] = i;
            t1[1] = i + 1;
            rel.insert(t1);
            rel.insert(t1);

            t1[0] = i;
            t1[1] = i;
            rel.insert(t1);
            rel.insert(t1);
        }

    for (uint64_t ep = 0; ep < epochs; ++ep)
        for (uint64_t j = 0; j < max; ++j)
        {
            const uint64_t i = ep*0x300000 + j;

            t1[0] = i;
            t1[1] = i + 1;
            rel.insert(t1);
            rel.insert(t1);

            t1[0] = i;
            t1[1] = i;
            rel.insert(t1);
            rel.insert(t1);
        }

    int rel_count = 0;
    t1[0] = -1; t1[1] = -1;
    tuple<2> selectall(t1);
    for (relation<2>::iter it(rel, selectall); it.more(); it.advance())
      rel_count++;

    assert(rel_count == epochs * max * 2);
    */

}

// Driver program to test above functions
int main()
{
    u64 start = utime();
    testseq();
    std::cout << ((utime() - start)/1000000) << std::endl;

    return 0;
}
