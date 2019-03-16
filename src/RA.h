#ifndef __RA__
#define __RA__

#include <iostream>
#include "relation.h"
#include "compat.h"


//template<unsigned arity>
static relation<2> join(relation<2>& r1, relation<2>& r2, int lc);

#if 1
static bool union1(relation<2>& r1, relation<2>& r2, int lc);

static bool union1(relation<2>& r1, relation<2>& r2, int lc)
{
    int count = 0;
    tuple<2> t;
    t[0] = -1;
    t[1] = -1;
    tuple<2> selectall(t);

    //std::cout << std::endl;
    //std::cout << "K" << std::endl;
    for (relation<2>::iter it(r2, selectall); it.more(); it.advance())
    {
      //std::cout << (*it)[0];
      //for (unsigned i = 1; i < 2; ++i)
      //    std::cout << ", " << (*it)[i];
      //std::cout << std::endl;

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

static relation<2> join(relation<2>& delT, relation<2>& G, int lc)
{
    tuple<2> t;
    t[0] = -1;
    t[1] = -1;
    tuple<2> selectall(t);

    relation<2> delTT;

    //std::cout << "T relation" << std::endl;
    for (relation<2>::iter dit(delT, selectall); dit.more(); dit.advance())
    {
      //std::cout << (*it)[0];
      //for (unsigned i = 1; i < 2; ++i)
      //    std::cout << ", " << (*it)[i];
      //std::cout << std::endl;

      tuple<2> s;
      s[0] = (*dit)[1];
      s[1] = -1;
      tuple<2> select(s);

      for (relation<2>::iter git(G, select); git.more(); git.advance())
      {
        //std::cout << "D" << (*git)[0];
        //for (unsigned i = 1; i < 2; ++i)
        //    std::cout << ", " << (*git)[i];
        //std::cout << std::endl;

        tuple<2> dt;
        dt[0] = (*dit)[0];
        dt[1] = (*git)[1];

        //if (lc == 1)
        //  std::cout << dt[0] << "\t" << dt[1] << std::endl;

        delTT.insert(dt);
      }
    }

    //if (lc == 0)
    //for (relation<2>::iter Tit(delTT, selectall); Tit.more(); Tit.advance())
    //std::cout << (*Tit)[0] << "\t" << (*Tit)[1] << "\tA" << std::endl;

    return delTT;
}



#endif
