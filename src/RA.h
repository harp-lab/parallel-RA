#ifndef __RA__
#define __RA__

#include <iostream>
#include "btree_relation.h"
#include "compat.h"


uint64_t utime()
{
    return std::chrono::high_resolution_clock::now().time_since_epoch().count();
}

relation<2> * join(relation<2>* r1, relation<2>& r2, relation<2>& r3, int lc, int* lb, int* rtc, u64* time);

inline static u64 tunedhash(const u8* bp, const u32 len)
{
    u64 h0 = 0xb97a19cb491c291d;
    u64 h1 = 0xc18292e6c9371a17;
    const u8* const ep = bp+len;
    while (bp < ep)
    {
        h1 ^= *bp;
        h1 *= 31;
        h0 ^= (((u64)*bp) << 17) ^ *bp;
        h0 *= 0x100000001b3;
        h0 = (h0 >> 7) | (h0 << 57);
        ++bp;
    }

    return h0 ^ h1;// ^ (h1 << 31);
}



static u64 outer_hash(const u64 val)
{
    return tunedhash((u8*)(&val),sizeof(u64));
}

#if 0
static void file_io(const char* filename, relation<2>& T);


static bool union1(relation<2>& r1, relation<2>& r2, int lc);


static bool union1(relation<2>& r1, relation<2>& r2, int lc)
{
    int count = 0;
    tuple<2> t;
    t[0] = -1;
    t[1] = -1;
    tuple<2> selectall(t);


    for (relation<2>::iter it(r2, selectall); it.more(); it.advance())
    {
      tuple<2> t1;
      t1[0] = (*it)[0];
      t1[1] = (*it)[1];

      if (r1.insert(t1) == true)
      {
        if (t1[0] == 2 && t1[1] == 22)
        {
          std::cout << lc << " Shout T" << std::endl;

          //for (relation<2>::iter it2(r1, selectall); it2.more(); it2.advance())
          //  std::cout << (*it2)[1] << "\t" << (*it2)[0] << std::endl;
        }
        count++;
      }
      else
      {
        if (t1[0] == 2 && t1[1] == 22)
          std::cout << lc << " Shout F" << std::endl;
      }
    }

    if (count == 0)
      return false;
    else
      return true;
}
#endif




/*
void file_io(const char* filename, relation<2>& T)
{
    tuple<2> t;
    t[0] = -1; t[1] = -1;
    tuple<2> selectall(t);

    std::ofstream myfile;
    myfile.open (filename);
    for (relation<2>::iter Tit(T, selectall); Tit.more(); Tit.advance())
      myfile << (*Tit)[0] << "\t" << (*Tit)[1] << "\n";
    myfile.close();
}
*/

#endif
