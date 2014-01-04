//#####################################################################
// Copyright 2009, Eftychios Sifakis, Yongning Zhu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RANGE_ITERATOR
//#####################################################################
#ifndef __RANGE_ITERATOR__
#define __RANGE_ITERATOR__

namespace PhysBAM{

template<int d,int stride=1> class RANGE_ITERATOR;

template<int d> 
class RANGE_ITERATOR<d,1>
{
    typedef VECTOR<int,d> TV_INT;
    typedef RANGE<TV_INT> T_RANGE;

    const T_RANGE& range;
    TV_INT index;

public:
    RANGE_ITERATOR(const T_RANGE& range_input)
        :range(range_input)
    {
        Reset();
    }

    void Reset()
    {index=range.min_corner;}

    bool Valid() const
    {return index.x<=range.max_corner.x;}

    void Next()
    {for(int i=d;i>=1;i--) if(index(i)<range.max_corner(i) || i==1){index(i)++;return;} else index(i)=range.min_corner(i);}

    const TV_INT& Index()
    {return index;}
};

template<int d,int stride>
class RANGE_ITERATOR
{
    STATIC_ASSERT((stride!=1));
    typedef VECTOR<int,d> TV_INT;
    typedef RANGE<TV_INT> T_RANGE;

    const T_RANGE& range;
    TV_INT index;

public:
    RANGE_ITERATOR(const T_RANGE& range_input)
        :range(range_input)
    {
        Reset();
    }

    void Reset()
    {index=range.min_corner;}

    bool Valid() const
    {return index.x<=range.max_corner.x;}

    void Next()
    {for(int i=d;i>=1;i--) if(index(i)+stride<=range.max_corner(i) || i==1){index(i)+=stride;return;} else index(i)=range.min_corner(i);}

    const TV_INT& Index()
    {return index;}
};

//#####################################################################
}
#endif
