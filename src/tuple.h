#ifndef __TUPLE__
#define __TUPLE__

template<unsigned arity>
class tuple
{
    typedef tuple<arity-1> subtuple;
private:
    s64 vals[arity];

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
    }

    s64& operator[](unsigned i)
    {
        return vals[i];
    }

    const s64& operator[](unsigned i) const
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


#endif
