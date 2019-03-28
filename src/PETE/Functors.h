// -*- C++ -*-
// ACL:license
// ----------------------------------------------------------------------
// This software and ancillary information (herein called "SOFTWARE")
// called PETE (Portable Expression Template Engine) is
// made available under the terms described here.  The SOFTWARE has been
// approved for release with associated LA-CC Number LA-CC-99-5.
//
// Unless otherwise indicated, this SOFTWARE has been authored by an
// employee or employees of the University of California, operator of the
// Los Alamos National Laboratory under Contract No.  W-7405-ENG-36 with
// the U.S. Department of Energy.  The U.S. Government has rights to use,
// reproduce, and distribute this SOFTWARE. The public may copy, distribute,
// prepare derivative works and publicly display this SOFTWARE without
// charge, provided that this Notice and any statement of authorship are
// reproduced on all copies.  Neither the Government nor the University
// makes any warranty, express or implied, or assumes any liability or
// responsibility for the use of this SOFTWARE.
//
// If SOFTWARE is modified to produce derivative works, such modified
// SOFTWARE should be clearly marked, so as not to confuse it with the
// version available from LANL.
//
// For more information about PETE, send e-mail to pete@acl.lanl.gov,
// or visit the PETE web page at http://www.acl.lanl.gov/pete/.
// ----------------------------------------------------------------------
// ACL:license

#ifndef PETE_PETE_FUNCTORS_H
#define PETE_PETE_FUNCTORS_H

///////////////////////////////////////////////////////////////////////////////
//
// WARNING: THIS FILE IS FOR INTERNAL PETE USE. DON'T INCLUDE IT YOURSELF
//
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Include files.
//-----------------------------------------------------------------------------

#include <iterator>

namespace qmcplusplus
{
//-----------------------------------------------------------------------------
//
// CLASS NAME
//   LeafFunctor<LeafType, LeafTag>
//
// DESCRIPTION
//   LeafType is the type of something at the leaf of the expression tree.
//   LeafTag specifies the operation being applied.
//
//   LeafFunctors are used by ForEach to apply operations to the leaves of the
//   expression tree. Typical functors are evaluators, counters, etc.
//   Users define functors for use with ForEach by specializing
//   the struct LeafFunctor<LeafType, LeafTag> for the user defined Functor and
//   any Leaf types that are necessary.
//
//   This isn't a functor in the conventional sense since it isn't invoked
//   operator() on a LeafFunctor object. Instead, the static apply member
//   function is called without an object. In a lot of ways, this is better
//   and more flexible than a regular functor.
//
//   LeafFunctor specializations must define the following:
//
//      typedef ... Type_t;
//         - the return type of the functor.
//      static Type_t apply(const LeafType &l, const LeafTag &f) {}
//         - evaluates the functor on leaf l.
//
//-----------------------------------------------------------------------------

template<class LeafType, class LeafTag>
struct LeafFunctor
  { };


//-----------------------------------------------------------------------------
//
// CLASS NAMES
//   EvalLeaf1-7, LeafFunctor<Scalar<T>, EvalLeaf1-7 >
//
// DESCRIPTION
//   EvalLeaf1 through EvalLeaf7 are used to evaluate leaves using 1 to 7
//   integer indices. We supply the tags and scalar versions here. Users must
//   supply specializations for their own containers.
//
//-----------------------------------------------------------------------------

// 1D

struct EvalLeaf1
{
  int i1_m;
  inline EvalLeaf1(int i1) : i1_m(i1) { }
  inline int val1() const
  {
    return i1_m;
  }
};

template<class T>
struct LeafFunctor<Scalar<T>, EvalLeaf1>
{
  typedef T Type_t;
  inline static
  const Type_t &apply(const Scalar<T> &s, const EvalLeaf1 &)
  {
    return s.value();
  }
};

// 2D

struct EvalLeaf2
{
  int i1_m, i2_m;
  inline EvalLeaf2(int i1, int i2) : i1_m(i1), i2_m(i2) { }
  inline int val1() const
  {
    return i1_m;
  }
  inline int val2() const
  {
    return i2_m;
  }
};

template<class T>
struct LeafFunctor<Scalar<T>, EvalLeaf2>
{
  typedef T Type_t;
  inline static
  const Type_t &apply(const Scalar<T> &s, const EvalLeaf2 &)
  {
    return s.value();
  }
};

// 3D

struct EvalLeaf3
{
  int i1_m, i2_m, i3_m;
  inline EvalLeaf3(int i1, int i2, int i3)
    : i1_m(i1), i2_m(i2), i3_m(i3) { }
  inline int val1() const
  {
    return i1_m;
  }
  inline int val2() const
  {
    return i2_m;
  }
  inline int val3() const
  {
    return i3_m;
  }
};

template<class T>
struct LeafFunctor<Scalar<T>, EvalLeaf3>
{
  typedef T Type_t;
  inline static
  const Type_t &apply(const Scalar<T> &s, const EvalLeaf3 &)
  {
    return s.value();
  }
};

// 4D

struct EvalLeaf4
{
  int i1_m, i2_m, i3_m, i4_m;
  inline EvalLeaf4(int i1, int i2, int i3, int i4)
    : i1_m(i1), i2_m(i2), i3_m(i3), i4_m(i4) { }
  inline int val1() const
  {
    return i1_m;
  }
  inline int val2() const
  {
    return i2_m;
  }
  inline int val3() const
  {
    return i3_m;
  }
  inline int val4() const
  {
    return i4_m;
  }
};

template<class T>
struct LeafFunctor<Scalar<T>, EvalLeaf4>
{
  typedef T Type_t;
  inline static
  const Type_t &apply(const Scalar<T> &s, const EvalLeaf4 &)
  {
    return s.value();
  }
};

// 5D

struct EvalLeaf5
{
  int i1_m, i2_m, i3_m, i4_m, i5_m;
  inline EvalLeaf5(int i1, int i2, int i3, int i4, int i5)
    : i1_m(i1), i2_m(i2), i3_m(i3), i4_m(i4), i5_m(i5) { }
  inline int val1() const
  {
    return i1_m;
  }
  inline int val2() const
  {
    return i2_m;
  }
  inline int val3() const
  {
    return i3_m;
  }
  inline int val4() const
  {
    return i4_m;
  }
  inline int val5() const
  {
    return i5_m;
  }
};

template<class T>
struct LeafFunctor<Scalar<T>, EvalLeaf5>
{
  typedef T Type_t;
  inline static
  const Type_t &apply(const Scalar<T> &s, const EvalLeaf5 &)
  {
    return s.value();
  }
};

// 6D

struct EvalLeaf6
{
  int i1_m, i2_m, i3_m, i4_m, i5_m, i6_m;
  inline EvalLeaf6(int i1, int i2, int i3, int i4, int i5, int i6)
    : i1_m(i1), i2_m(i2), i3_m(i3), i4_m(i4), i5_m(i5), i6_m(i6) { }
  inline int val1() const
  {
    return i1_m;
  }
  inline int val2() const
  {
    return i2_m;
  }
  inline int val3() const
  {
    return i3_m;
  }
  inline int val4() const
  {
    return i4_m;
  }
  inline int val5() const
  {
    return i5_m;
  }
  inline int val6() const
  {
    return i6_m;
  }
};

template<class T>
struct LeafFunctor<Scalar<T>, EvalLeaf6>
{
  typedef T Type_t;
  inline static
  const Type_t &apply(const Scalar<T> &s, const EvalLeaf6 &)
  {
    return s.value();
  }
};

// 7D

struct EvalLeaf7
{
  int i1_m, i2_m, i3_m, i4_m, i5_m, i6_m, i7_m;
  inline EvalLeaf7(int i1, int i2, int i3, int i4, int i5, int i6,
                   int i7)
    : i1_m(i1), i2_m(i2), i3_m(i3), i4_m(i4), i5_m(i5), i6_m(i6), i7_m(i7) { }
  inline int val1() const
  {
    return i1_m;
  }
  inline int val2() const
  {
    return i2_m;
  }
  inline int val3() const
  {
    return i3_m;
  }
  inline int val4() const
  {
    return i4_m;
  }
  inline int val5() const
  {
    return i5_m;
  }
  inline int val6() const
  {
    return i6_m;
  }
  inline int val7() const
  {
    return i7_m;
  }
};

template<class T>
struct LeafFunctor<Scalar<T>, EvalLeaf7>
{
  typedef T Type_t;
  inline static
  const Type_t &apply(const Scalar<T> &s, const EvalLeaf7 &)
  {
    return s.value();
  }
};


//-----------------------------------------------------------------------------
//
// CLASS NAME
//   IncrementLeaf, LeafFunctor<{T, Scalar<T>}, IncrementLeaf >
//
// DESCRIPTION
//   A leaf-tag and functor used to increment an iterator for scalars and
//   leaves made up of STL iterators.
//
//-----------------------------------------------------------------------------

struct IncrementLeaf
  { };

template<class T>
struct LeafFunctor<T, IncrementLeaf>
{
  typedef int Type_t;
  inline static
  Type_t apply(const T &cl, const IncrementLeaf &)
  {
    T &l = const_cast<T &>(cl);
    ++l;
    return 0;
  }
};

#if defined(__MWERKS__)

// Workaround for screwy CWPro 4.1 bug.

template <class T>
struct LeafFunctor<const T*, IncrementLeaf>
{
  typedef int Type_t;
  inline static
  Type_t apply(const T* & const ci, const IncrementLeaf &)
  {
    T* &i = const_cast<T* &>(ci);
    ++i;
    return 0;
  }
};

#endif

template<class T>
struct LeafFunctor<Scalar<T>, IncrementLeaf>
{
  typedef int Type_t;
  inline static
  Type_t apply(const Scalar<T> &, const IncrementLeaf &)
  {
    return 0;
  }
};


//-----------------------------------------------------------------------------
//
// CLASS NAME
//   DecrementLeaf, LeafFunctor<{T, Scalar<T>}, DecrementLeaf >
//
// DESCRIPTION
//   A leaf-tag and functor used to decrement an iterator for scalars and
//   leaves made up of STL iterators.
//
//-----------------------------------------------------------------------------

struct DecrementLeaf
  { };

template<class T>
struct LeafFunctor<T, DecrementLeaf>
{
  typedef int Type_t;
  inline static
  Type_t apply(const T &cl, const DecrementLeaf &)
  {
    T &l = const_cast<T &>(cl);
    --l;
    return 0;
  }
};

#if defined(__MWERKS__)
// Workaround for screwy CWPro 4.1 bug.
template <class T>
struct LeafFunctor<const T*, DecrementLeaf>
{
  typedef int Type_t;
  inline static
  Type_t apply(const T* & const ci, const IncrementLeaf &)
  {
    T* &i = const_cast<T* &>(ci);
    --i;
    return 0;
  }
};
#endif

template<class T>
struct LeafFunctor<Scalar<T>, DecrementLeaf>
{
  typedef int Type_t;
  inline static
  Type_t apply(const Scalar<T> &, const DecrementLeaf &)
  {
    return 0;
  }
};

//-----------------------------------------------------------------------------
//
// CLASS NAME
//   DereferenceLeaf, LeafFunctor<{T, Scalar<T>}, DereferenceLeaf >
//
// DESCRIPTION
//   A leaf-tag and functor used to dereference an iterator for scalars and
//   leaves made up of STL iterators.
//
//-----------------------------------------------------------------------------

struct DereferenceLeaf
  { };

template<class ForwardIterator>
struct LeafFunctor<ForwardIterator, DereferenceLeaf>
{
  typedef typename std::iterator_traits<ForwardIterator>::value_type Type_t;
  inline static
  Type_t apply(const ForwardIterator &i, const DereferenceLeaf &)
  {
    return *i;
  }
};

#if defined(__MWERKS__)
// Workaround for screwy CWPro 4.1 bug.
template <class T>
struct LeafFunctor<const T*, DereferenceLeaf>
{
  typedef T Type_t;
  inline static
  Type_t apply(const T *i, const DereferenceLeaf &)
  {
    return *i;
  }
};
#endif

template<class T>
struct LeafFunctor<Scalar<T>, DereferenceLeaf>
{
  typedef T Type_t;
  inline static
  const Type_t &apply(const Scalar<T> &s, const DereferenceLeaf &)
  {
    return s.value();
  }
};

}

#endif // PETE_PETE_FUNCTORS_H

// ACL:rcsinfo
