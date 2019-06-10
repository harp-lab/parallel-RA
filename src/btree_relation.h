#ifndef __BTREE_RELATION__
#define __BTREE_RELATION__

#include <iostream>
#include "tuple.h"
#include "btree.h"
#include "compat.h"


template<unsigned arity>
class relation
{
    typedef relation<arity-1> subrelation;
    typedef BTree<subrelation> btree;

private:
    btree tree;

public:
    relation()
    {
    }

    relation(const relation<arity>& r)
    {
        (*this).tree = r.tree;
    }

    ~relation()
    {
        for (typename btree::iter it(&tree, -1); it.more(); it.advance())
        {
            subrelation *sr = it.getval();
            if (sr != (void*)1)
                delete sr;
        }
    }


    class iter
    {
    private:
        tuple<arity> selection;
        typename btree::iter btit;
        typename subrelation::iter srit;

    public:

        iter(const relation<arity>& r, const tuple<arity>& sel)
            : selection(sel),
              btit(&r.tree, selection[0]),
              srit(btit.getval()->select(selection.tail()))
        {
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
                    btit.advance();
                    if (!btit.more())
                      return;

                    srit = btit.getval()->select(selection.tail());
                }
            }
        }

        bool more()
        {
            return srit.more();
        }

        tuple<arity> operator *()
        {
            return tuple<arity>(btit.getkey(), *srit);
        }

        void operator =(relation<arity>::iter it)
        {
            selection = it.selection;
            btit = it.btit;
            srit = it.srit;
        }
    };

    void operator =(relation<arity> r)
    {
        tree = r.tree;
    }

    relation(u32 buffer_size, u64* buffer)
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
        subrelation* sr = tree.search(t[0]);
        bool modified = false;
        if (sr == NULL)
        {
            modified = true;
            sr = new subrelation();
            tree.insert(t[0], sr);
        }

        return sr->insert(t.tail()) || modified;
    }


    iter select(const tuple<arity>& sel)
    {
        iter it(*this, sel);
        return it;
    }
};


template<>
class relation<1>
{
    typedef BTree<void> btree;

private:
    btree tree;

public:

    relation()
    {
    }

    relation(const relation<1>& r)
    {
        (*this).tree = r.tree;
    }

    class iter
    {
    private:
        tuple<1> selection;
        typename btree::iter btit;

    public:
        iter(const relation<1>& r, const tuple<1>& sel)
            : selection(sel), btit(&r.tree, selection[0])
        {
        }

        void advance()
        {
            btit.advance();
        }

        bool more()
        {
            return btit.more();
        }

        tuple<1> operator *()
        {
            tuple<1> t;
            t[0] = btit.getkey();
            return t;
        }

        void operator =(relation<1>::iter it)
        {
            selection = it.selection;
            btit = it.btit;
        }
    };

    bool insert(const tuple<1> t)
    {
        return tree.insert(t[0], (void*)1);
    }

    iter select(const tuple<1>& sel)
    {
        iter it(*this, sel);
        return it;
    }
};

#endif
