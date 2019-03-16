class tuple
{
private:
    u64 vals[2];

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

    tuple(u64 n, const tuple<arity-1>& t)
    {
        vals[0] = n;
        for (int i = 1; i < arity; i++)
            vals[i] = t.vals[i];
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
            vals[i] = t.vals[i];
    }

    tuple<arity-1> tail()
    {
        tuple<arity-1> st;
        for (unsigned i = 1; i < arity; ++i)
            st[i-1] = operator[](i);
        return st;
    }
};


template<unsigned arity>
class tuple
{
private:
    u64 vals[2];

public:
    tuple()
    {
        for (u32 i = 0; i < 1; i++)
            vals[i] = -1;
    }

    tuple(const tuple<1>& t)
    {
        (*this) = t;
    }

    tuple(u64 n, const tuple<0>& t)
    {
        vals[0] = n;
        for (int i = 1; i < 1; i++)
            vals[i] = t.vals[i];
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
        for (u32 i = 0; i < 1; i++)
            vals[i] = t.vals[i];
    }

    tuple<0> tail()
    {
        tuple<0> st;
        return st;
    }
};
