#ifndef __RELATION__
#define __RELATION__

#include <iostream>
#include "btree.h"
#include "compat.h"


template<unsigned arity>
class tuple
{
    typedef tuple<arity-1> subtuple;
private:
    u64 vals[arity];

public:
    tuple()
    {
        for (u32 i = 0; i < arity; i++)
            vals[i] = -1;
    }

    tuple(const tuple<arity>& t)
    {
        (*this) = t;
    }

    tuple(u64 n, const subtuple& t)
    {
        vals[0] = n;
        for (u32 i = 1; i < arity; i++)
            vals[i] = t[i - 1];

        //std::cout << "[2] key" << t[0] << std::endl;
    }

    u64& operator[](unsigned i)
    {
        return vals[i];
    }

    const u64& operator[](unsigned i) const
    {
        return vals[i];
    }

    void operator=(const tuple<arity>& t)
    {
        for (u32 i = 0; i < arity; i++)
            vals[i] = t[i];
    }


    subtuple tail() const
    {
        subtuple st;
        for (u32 i = 1; i < arity; ++i)
            st[i-1] = operator[](i);
        return st;
    }
};

/*
template<>
class tuple<1>
{
    //typedef tuple<0> subtuple;
private:
    u64 vals[1];

public:
    tuple()
    {
        vals[0] = -1;
    }

    tuple(const tuple<1>& t)
    {
        (*this) = t;
    }

    u64& operator[](unsigned i)
    {
        return vals[i];
    }

    const u64& operator[](unsigned i) const
    {
        return vals[i];
    }

    void operator=(const tuple<1>& t)
    {
        vals[0] = t[0];
    }


    tuple<0> tail()
    {
        tuple<0> st;
        //for (unsigned i = 1; i < arity; ++i)
        //    st[i-1] = operator[](i);
        return st;
    }
};
*/


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
        (*this) = r;
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
#ifdef LOGGING
            std::cout<<"[Constructor 2]" << btit.getkey() << std::endl;
#endif
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
#ifdef LOGGING
                    std::cout << "XXXX MORE" << btit.getkey() << std::endl;
#endif
                    btit.advance();
#ifdef LOGGING
                    std::cout << "YYYY MORE " << btit.getkey() << std::endl;
#endif
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
            //std::cout << "[1] key" << btit.getkey() << std::endl;
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
#ifdef LOGGING
        std::cout << "[Relation] [" << arity << "] " << "Search: " << t[0] << std::endl;
#endif
        subrelation* sr = tree.search(t[0]);
        bool modified = false;
        if (sr == NULL)
        {
#ifdef LOGGING
            std::cout << "[Relation] [" << arity << "] " << "Insert: " << t[0] << std::endl;
#endif
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

    class iter
    {
    private:
        tuple<1> selection;
        typename btree::iter btit;

    public:
        iter(const relation<1>& r, const tuple<1>& sel)
            : selection(sel), btit(&r.tree, selection[0])
        {
#ifdef LOGGING
            std::cout << "[Constructor 1] [iter 1] init [1]" << btit.getkey()  <<std::endl;
#endif
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
#ifdef LOGGING
        std::cout << "[Relation] [1] insert: " << t[0] << std::endl;
#endif
        return tree.insert(t[0], (void*)1);
    }

    iter select(const tuple<1>& sel)
    {
        iter it(*this, sel);
        return it;
    }
};

#endif
