#ifndef __RA__
#define __RA__

#include <iostream>
#include "relation.h"
#include "compat.h"

uint64_t utime()
{
    return std::chrono::high_resolution_clock::now().time_since_epoch().count();
}

relation<2> join(relation<2>& r1, relation<2>& r2, relation<2>& r3, int lc, int* lb, int* rtc, u64* time);

std::vector<u64*> joinV(std::vector<u64*>& delT, relation<2>& G, relation<2>& T, int lc, int* lb, int* running_t_count, u64* running_time);

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

relation<2> join(relation<2>& delT, relation<2>& G, relation<2>& T, int lc, int* lb, int* running_t_count, u64* running_time)
{
    u64 start = utime();
    tuple<2> t;
    t[0] = -1;
    t[1] = -1;
    tuple<2> selectall(t);

    relation<2> delTT;

    int count = 0;
    int tcount = 0;

    for (relation<2>::iter dit(delT, selectall); dit.more(); dit.advance())
    {
      tuple<2> s;
      s[0] = (*dit)[1];
      s[1] = -1;
      tuple<2> select(s);

      for (relation<2>::iter git(G, select); git.more(); git.advance())
      {
        tuple<2> dt;
        dt[0] = (*dit)[0];
        dt[1] = (*git)[1];

        if (T.insert(dt, lc) == true)
        {
            tcount++;
            if (delTT.insert(dt, lc) == true)
              count++;
        }
      }
    }

    if (count == 0)
        *lb = 1;
    else
        *lb = 0;

    *running_t_count = *running_t_count + tcount;

    u64 end = utime();
    u64 dTime = (end - start) / 1000000;
    *running_time = *running_time + dTime;

    std::cout << lc << " [" << dTime << "]  [" << *running_time << "] : Delta count " << count << " T count: " << *running_t_count << " : " << std::endl;

    return delTT;
}


std::vector<u64> joinV(std::vector<u64>& delT, relation<2>& G, relation<2>& T, int lc, int* lb, int* running_t_count, u64* running_time)
{
    //std::vector<u64> r;

    u64 start = utime();
    tuple<2> t;
    t[0] = -1;
    t[1] = -1;
    tuple<2> selectall(t);

    std::vector<u64> delTT;

    int count = 0;
    int tcount = 0;

    for (std::vector<u64>::iterator it = delT.begin() ; it != delT.end(); ++it)
    {
      u64 a = *it;
      ++it;
      u64 b = *it;

      tuple<2> s;
      s[0] = b;
      s[1] = -1;
      tuple<2> select(s);

      for (relation<2>::iter git(G, select); git.more(); git.advance())
      {
        tuple<2> dt;
        dt[0] = a;
        dt[1] = (*git)[1];

        if (T.insert(dt, lc) == true)
        {
            tcount++;
            delTT.push_back(a);
            delTT.push_back(dt[1]);
            count++;
        }
      }
    }

    if (count == 0)
        *lb = 1;
    else
        *lb = 0;

    *running_t_count = *running_t_count + tcount;

    u64 end = utime();
    u64 dTime = (end - start) / 1000000;
    *running_time = *running_time + dTime;

    std::cout << lc << " [" << dTime << "]  [" << *running_time << "] : Delta count " << count << " T count: " << *running_t_count << " : " << std::endl;

    return delTT;
}

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
