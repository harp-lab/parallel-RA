#ifndef __GOOGLE_RELATION__
#define __GOOGLE_RELATION__


#include "compat.h"
#include "tuple.h"
#include "btree/btree_map.h"



//typedef btree::btree_map<u64, u64> Map0;
//typedef btree::btree_map<u64, Map0* > Map1;

class google_relation
{
    typedef google_relation subrelation;
    typedef btree::btree_map<u64, subrelation*> btree;

private:
    u32 arity;
    btree tree;

public:

    google_relation()
    {
    }

    void setArity(u32 a) {arity = a;}

    bool insert(const u64* t, u32 arity)
    {
        auto temp = tree.find(t[0]);
        bool modified = false;
        subrelation *sr;
        if (temp == tree.end())
        {
            modified = true;
            sr = new subrelation();
            tree.insert(std::make_pair(t[0], sr));
        }
        else
            sr = temp->second;

        return sr->insert(t++, arity - 1) || modified;
    }

    class iter
    {
    private:
        btree root_tree;
        const s64* selection;
        typename btree::const_iterator btit;
        typename subrelation::iter srit;


    public:

        iter(const google_relation& r, const s64* sel)
            : root_tree(r.tree),
              selection(sel),
              srit(btit->second->select(selection++))
        {
            if (selection[0] == -1)
                btit = r.tree.begin();
            else
                btit = (r.tree).find(selection[0]);

        }

        iter (const iter& it)
            : selection(it.selection),
              btit(it.btit),
              srit(it.srit)
        {
        }


        void advance()
        {
            if (srit.more())
            {
                srit.advance();
                while (!srit.more())
                {
                    btit++;
                    if (btit == root_tree.end())
                        return;

                    srit = btit->second->select(selection++);
                }
            }
        }

        bool more()
        {
            return srit.more();
        }


        u64* operator *()
        {
            return 0;
            //return tuple(btit->first, *srit);
        }


        void operator =(google_relation::iter it)
        {
        }


    };

    iter select(const s64* sel)
    {
        iter it(*this, sel);
        return it;
    }

    void operator =(google_relation r)
    {
        tree = r.tree;
    }

};

#if 0
template<unsigned arity>
class google_relation
{
    typedef google_relation<arity-1> subrelation;
    typedef btree::btree_map<u64, subrelation* > btree;

private:
    btree tree;

public:
    google_relation()
    {
    }

    google_relation(const google_relation<arity>& r)
    {
        (*this).tree = r.tree;
    }

    ~google_relation()
    {
        for(typename btree::iterator iy2 = tree.begin(); iy2 != tree.end(); iy2++)
            if (iy2->second != (void*)1)
                delete iy2->second;
    }

    //
    class iter
    {
    private:
        btree root_tree;
        tuple<arity> selection;
        typename btree::const_iterator btit;
        typename subrelation::iter srit;

    public:

        iter(const google_relation<arity>& r, const tuple<arity>& sel)
            : root_tree(r.tree),
              selection(sel),
              srit(btit->second->select(selection.tail()))
        {
            if (selection[0] == -1)
                btit = r.tree.begin();
            else
                btit = (r.tree).find(selection[0]);
        }

        iter (const iter& it)
            : selection(it.selection),
              btit(it.btit),
              srit(it.srit)
        {
        }

        void advance()
        {
            if (srit.more())
            {
                srit.advance();
                while (!srit.more())
                {
                    btit++;
                    if (btit == root_tree.end())
                        return;

                    srit = btit->second->select(selection.tail());
                }
            }
        }

        bool more()
        {
            return srit.more();
        }

        tuple<arity> operator *()
        {
            return tuple<arity>(btit->first, *srit);
        }

        void operator =(google_relation<arity>::iter it)
        {
            selection = it.selection;
            btit = it.btit;
            srit = it.srit;
        }
    };
    //

    void operator =(google_relation<arity> r)
    {
        tree = r.tree;
    }

    google_relation(u32 buffer_size, u64* buffer)
    {
        for (u32 i = 0; i < buffer_size; i = i + arity)
        {
            tuple<arity> t;
            for (unsigned j = 0; j < arity; j++)
                t[j] = buffer[i + j];
            insert(t);
        }
    }

    void initialize(u32 buffer_size, u64* buffer)
    {
        for (u32 i = 0; i < buffer_size; i = i + arity)
        {
            tuple<arity> t;
            for (unsigned j = 0; j < arity; j++)
                t[j] = buffer[i + j];
            insert(t);
        }
    }


    bool insert(const tuple<arity>& t)
    {
        auto temp = tree.find(t[0]);
        bool modified = false;
        subrelation *sr;
        if (temp == tree.end())
        {
            modified = true;
            sr = new subrelation();
            tree.insert(std::make_pair(t[0], sr));
        }
        else
            sr = temp->second;

        return sr->insert(t.tail()) || modified;
    }



    iter select(const tuple<arity>& sel)
    {
        iter it(*this, sel);
        return it;
    }
};


template<>
class google_relation<1>
{
    typedef btree::btree_map<u64, u64> btree;

private:
    btree tree;

public:

    google_relation()
    {
    }

    google_relation(const google_relation<1>& r)
    {
        (*this).tree = r.tree;
    }


    class iter
    {
    private:
        btree root_tree;
        tuple<1> selection;
        typename btree::const_iterator btit;

    public:
        iter(const google_relation<1>& r, const tuple<1>& sel)
            : root_tree(r.tree),
              selection(sel)
        {
            if (selection[0] == -1)
                btit = r.tree.begin();
            else
                btit = r.tree.find(selection[0]);
        }

        void advance()
        {
            btit++;
        }

        bool more()
        {
            return btit == root_tree.end();
        }

        tuple<1> operator *()
        {
            tuple<1> t;
            t[0] = btit->first;
            return t;
        }

        void operator =(google_relation<1>::iter it)
        {
            selection = it.selection;
            btit = it.btit;
        }
    };


    bool insert(const tuple<1> t)
    {
        auto it = tree.find(t[0]);
        if ( it != tree.end() ) {
            return false;
        }
        else
        {
            tree.insert(std::make_pair(t[0], 0));
            return true;
        }
    }

    iter select(const tuple<1>& sel)
    {
        iter it(*this, sel);
        return it;
    }
};
#endif

#if 0
typedef btree::btree_map<u64, u64> Map0;
typedef btree::btree_map<u64, btree::btree_map<u64, u64>* > Map1;

class google_relation
{

private:
    Map1 tree;

public:
    class iter
    {
    private:
        tuple<2> selection;
        typename Map1::iterator M1_it;
        typename Map0::iterator M0_it;

    public:
    };

    google_relation()
    {

    }

    google_relation(u32 buffer_size, u64* buffer)
    {
        for (u32 i = 0; i < buffer_size; i = i + 2)
        {
            tuple<2> t;
            for (unsigned j = 0; j < 2; j++)
                t[j] = buffer[i + j];
            insert(t);
        }
    }


    bool insert(const tuple<2>& t)
    {
        auto it = tree.find(t[0]);
        if( it != tree.end())
        {
            auto it2 = (it->second)->find(t[1]);
            if( it2 == (it->second)->end() )
            {
                (it->second)->insert(std::make_pair(t[1], 0));
                tree[t[0]] = it->second;
                return true;
            }
        }
        else
        {
            Map0 *k = new Map0;
            k->insert(std::make_pair(t[1], 0));
            tree.insert(std::make_pair(t[0],k));
            return true;
        }

        return false;
    }


    bool insert(const u64& t1, const u64& t2)
    {
        auto it = tree.find(t1);
        if( it != tree.end())
        {
            auto it2 = (it->second)->find(t2);
            if( it2 == (it->second)->end() )
            {
                (it->second)->insert(std::make_pair(t2, 0));
                tree[t1] = it->second;
                return true;
            }
        }
        else
        {
            Map0 *k = new Map0;
            k->insert(std::make_pair(t2, 0));
            tree.insert(std::make_pair(t1,k));
            return true;
        }

        return false;
    }


    //iter select(const tuple<2>& sel)
    //{

    //}
};
#endif

#endif
