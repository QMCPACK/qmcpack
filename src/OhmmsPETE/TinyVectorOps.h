//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign   
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#ifndef OHMMS_TINYVECTOR_OPERATORS_H
#define OHMMS_TINYVECTOR_OPERATORS_H
#include <complex>

namespace qmcplusplus
{

template<class T1>
struct BinaryReturn<T1, std::complex<T1>, OpMultiply >
{
  typedef std::complex<T1> Type_t;
};

template<class T1>
struct BinaryReturn<std::complex<T1>, T1, OpMultiply >
{
  typedef std::complex<T1> Type_t;
};

///////////////////////////////////////////////////////////////////////
//
// Assignment operators
// template<class T1, class T2, class OP> struct OTAssign {};
//
///////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Specializations for TinyVectors of arbitrary size.
//////////////////////////////////////////////////////////////////////
template<class T1, class T2, class OP, unsigned D>
struct OTAssign< TinyVector<T1,D> , TinyVector<T2,D> , OP >
{
  inline static void
  apply( TinyVector<T1,D>& lhs, const TinyVector<T2,D>& rhs, OP op)
  {
    for (unsigned d=0; d<D; ++d)
      op(lhs[d] , rhs[d]);
  }
};

template<class T1, class T2, class OP, unsigned D>
struct OTAssign< TinyVector<T1,D> , T2 , OP >
{
  inline static void
  apply( TinyVector<T1,D>& lhs, const T2& rhs, OP op )
  {
    for (unsigned d=0; d<D; ++d)
      op(lhs[d] , rhs);
  }
};

//////////////////////////////////////////////////////////////////////
// Specializations for TinyVectors with D=1.
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct OTAssign< TinyVector<T1,1> , TinyVector<T2,1> , OP >
{
  inline static void
  apply( TinyVector<T1,1>& lhs, const TinyVector<T2,1>& rhs, OP op)
  {
    op(lhs[0] , rhs[0] );
  }
};

template<class T1, class T2, class OP>
struct OTAssign< TinyVector<T1,1> , T2 , OP >
{
  inline static void
  apply( TinyVector<T1,1>& lhs, const T2& rhs, OP op )
  {
    op(lhs[0] , rhs );
  }
};

//////////////////////////////////////////////////////////////////////
// Specializations for TinyVectors with D=2.
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct OTAssign< TinyVector<T1,2> , TinyVector<T2,2> , OP >
{
  inline static void
  apply( TinyVector<T1,2>& lhs, const TinyVector<T2,2>& rhs, OP op)
  {
    op(lhs[0] , rhs[0] );
    op(lhs[1] , rhs[1] );
  }
};

template<class T1, class T2, class OP>
struct OTAssign< TinyVector<T1,2> , T2 , OP >
{
  inline static void
  apply( TinyVector<T1,2>& lhs, const T2& rhs, OP op )
  {
    op(lhs[0] , rhs);
    op(lhs[1] , rhs);
  }
};

//////////////////////////////////////////////////////////////////////
// Specializations for TinyVectors with D=3.
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct OTAssign< TinyVector<T1,3> , TinyVector<T2,3> , OP >
{
  inline static void
  apply( TinyVector<T1,3>& lhs, const TinyVector<T2,3>& rhs, OP op)
  {
    op(lhs[0] , rhs[0] );
    op(lhs[1] , rhs[1] );
    op(lhs[2] , rhs[2] );
  }
};

template<class T1, class T2, class OP>
struct OTAssign< TinyVector<T1,3> , TinyVector<std::complex<T2>,3> , OP >
{
  inline static void
  apply( TinyVector<T1,3>& lhs, const TinyVector<std::complex<T2>,3>& rhs, OP op)
  {
    op(lhs[0] , rhs[0].real() );
    op(lhs[1] , rhs[1].real() );
    op(lhs[2] , rhs[2].real() );
  }
};

template<class T1, class T2, class OP>
struct OTAssign< TinyVector<std::complex<T1>,3> , TinyVector<std::complex<T2>,3> , OP >
{
  inline static void
  apply( TinyVector<std::complex<T1>,3>& lhs, const TinyVector<std::complex<T2>,3>& rhs, OP op)
  {
    op(lhs[0] , rhs[0] );
    op(lhs[1] , rhs[1] );
    op(lhs[2] , rhs[2] );
  }
};

template<class T1, class T2, class OP>
struct OTAssign< TinyVector<T1,3> , T2 , OP >
{
  inline static void
  apply( TinyVector<T1,3>& lhs, const T2& rhs, OP op )
  {
    op(lhs[0] , rhs );
    op(lhs[1] , rhs );
    op(lhs[2] , rhs );
  }
};

///////////////////////////////////////////////////////////////////////
//
// Binary operators
//template<class T1, class T2, class OP> struct OTBinary {};
//
///////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Specializations for TinyVectors of arbitrary size.
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP, unsigned D>
struct OTBinary< TinyVector<T1,D> , TinyVector<T2,D> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyVector<Type_t,D>
  apply(const TinyVector<T1,D>& lhs, const TinyVector<T2,D>& rhs, OP op)
  {
    TinyVector<Type_t,D> ret;
    for (unsigned d=0; d<D; ++d)
      ret[d] = op(lhs[d] , rhs[d]);
    return ret;
  }
};

template<class T1, class T2, class OP, unsigned D>
struct OTBinary< TinyVector<T1,D> , T2 , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyVector<Type_t,D>
  apply(const TinyVector<T1,D>& lhs, const T2& rhs, OP op)
  {
    TinyVector<Type_t,D> ret;
    for (unsigned d=0 ; d<D; ++d)
      ret[d] = op(lhs[d] , rhs );
    return ret;
  }
};

template<class T1, class T2, class OP, unsigned D>
struct OTBinary< T1, TinyVector<T2,D> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyVector<Type_t,D>
  apply(const T1& lhs, const TinyVector<T2,D>& rhs, OP op)
  {
    TinyVector<Type_t,D> ret;
    for (unsigned d=0; d<D; ++d)
      ret[d] = op(lhs , rhs[d]);
    return ret;
  }
};

//////////////////////////////////////////////////////////////////////
// Specializations of OTBinary for TinyVectors with D=1.
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct OTBinary< TinyVector<T1,1> , TinyVector<T2,1> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyVector<Type_t,1>
  apply(const TinyVector<T1,1>& lhs, const TinyVector<T2,1>& rhs, OP op)
  {
    return TinyVector<Type_t,1>( op(lhs[0], rhs[0] ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< TinyVector<T1,1> , T2 , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyVector<Type_t,1>
  apply(const TinyVector<T1,1>& lhs, const T2& rhs, OP op)
  {
    return TinyVector<Type_t,1>( op(lhs[0], rhs ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< T1, TinyVector<T2,1> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyVector<Type_t,1>
  apply(const T1& lhs, const TinyVector<T2,1>& rhs, OP op)
  {
    return TinyVector<Type_t,1>( op(lhs, rhs[0] ) );
  }
};

//////////////////////////////////////////////////////////////////////
// Specializations of OTBinary for TinyVectors with D=2.
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct OTBinary< TinyVector<T1,2> , TinyVector<T2,2> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyVector<Type_t,2>
  apply(const TinyVector<T1,2>& lhs, const TinyVector<T2,2>& rhs, OP op)
  {
    return TinyVector<Type_t,2>( op(lhs[0], rhs[0]),op(lhs[1], rhs[1] ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< TinyVector<T1,2> , T2 , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyVector<Type_t,2>
  apply(const TinyVector<T1,2>& lhs, const T2& rhs, OP op)
  {
    return TinyVector<Type_t,2>( op(lhs[0], rhs ) ,op(lhs[1], rhs));
  }
};

template<class T1, class T2, class OP>
struct OTBinary< T1, TinyVector<T2,2> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyVector<Type_t,2>
  apply(const T1& lhs, const TinyVector<T2,2>& rhs, OP op)
  {
    return TinyVector<Type_t,2>( op(lhs, rhs[0] ) ,op(lhs, rhs[1] ) );
  }
};

//////////////////////////////////////////////////////////////////////
// Specializations of OTBinary for TinyVectors with D=3.
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct OTBinary< TinyVector<T1,3> , TinyVector<T2,3> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyVector<Type_t,3>
  apply(const TinyVector<T1,3>& lhs, const TinyVector<T2,3>& rhs, OP op)
  {
    return TinyVector<Type_t,3>( op(lhs[0], rhs[0] ) ,
                                 op(lhs[1], rhs[1] ) ,
                                 op(lhs[2], rhs[2] ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< TinyVector<T1,3> , T2 , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyVector<Type_t,3>
  apply(const TinyVector<T1,3>& lhs, const T2& rhs, OP op)
  {
    return TinyVector<Type_t,3>( op(lhs[0], rhs ) ,
                                 op(lhs[1], rhs ) ,
                                 op(lhs[2], rhs ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< T1, TinyVector<T2,3> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyVector<Type_t,3>
  apply(const T1& lhs, const TinyVector<T2,3>& rhs, OP op)
  {
    return TinyVector<Type_t,3>( op( lhs, rhs[0] ) ,
                                 op( lhs, rhs[1] ) ,
                                 op( lhs, rhs[2] ) );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations for TinyVector dot TinyVector
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, unsigned D>
struct OTDot< TinyVector<T1,D> , TinyVector<T2,D> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Type_t
  apply(const TinyVector<T1,D>& lhs, const TinyVector<T2,D>& rhs)
  {
    Type_t res = lhs[0]*rhs[0];
    for (unsigned d=1; d<D; ++d)
      res += lhs[d]*rhs[d];
    return res;
  }
};

template<class T1, class T2>
struct OTDot< TinyVector<T1,1> , TinyVector<T2,1> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Type_t
  apply(const TinyVector<T1,1>& lhs, const TinyVector<T2,1>& rhs)
  {
    return lhs[0]*rhs[0];
  }
};

template<class T1, class T2>
struct OTDot< TinyVector<T1,2> , TinyVector<T2,2> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Type_t
  apply(const TinyVector<T1,2>& lhs, const TinyVector<T2,2>& rhs)
  {
    return lhs[0]*rhs[0] + lhs[1]*rhs[1];
  }
};

template<class T1, class T2>
struct OTDot< TinyVector<T1,3> , TinyVector<T2,3> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Type_t
  apply(const TinyVector<T1,3>& lhs, const TinyVector<T2,3>& rhs)
  {
    return lhs[0]*rhs[0] + lhs[1]*rhs[1] + lhs[2]*rhs[2];
  }
};

template<class T1, class T2>
struct OTDot< TinyVector<T1,4> , TinyVector<T2,4> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Type_t
  apply(const TinyVector<T1,4>& lhs, const TinyVector<T2,4>& rhs)
  {
    return lhs[0]*rhs[0] + lhs[1]*rhs[1] + lhs[2]*rhs[2] + lhs[3]*rhs[3];
  }
};

/** specialization for real-complex TinyVector */
template<class T1>
struct OTDot< TinyVector<T1,3> , TinyVector<std::complex<T1>,3> >
{
  typedef T1 Type_t;
  inline static Type_t
  apply(const TinyVector<T1,3>& lhs, const TinyVector<std::complex<T1>,3>& rhs)
  {
    return lhs[0]*rhs[0].real() + lhs[1]*rhs[1].real() + lhs[2]*rhs[2].real();
  }
};

/** specialization for complex-real TinyVector */
template<class T1, class T2>
struct OTDot< TinyVector<std::complex<T1>,3> , TinyVector<T2,3> >
{
  typedef T1 Type_t;
  inline static Type_t
  apply(const TinyVector<std::complex<T1>,3>& lhs, const TinyVector<T2,3>& rhs)
  {
    return lhs[0].real()*rhs[0] + lhs[1].real()*rhs[1] + lhs[2].real()*rhs[2];
  }
};

/** specialization for complex-complex TinyVector */
template<class T1, class T2>
struct OTDot< TinyVector<std::complex<T1>,3> , TinyVector<std::complex<T2>,3> >
{
  typedef typename BinaryReturn<std::complex<T1>, std::complex<T2>, OpMultiply>::Type_t Type_t;
  inline static Type_t
  apply(const TinyVector<std::complex<T1>,3>& lhs, const TinyVector<std::complex<T2>,3>& rhs)
  {
    return std::complex<T1>(lhs[0].real()*rhs[0].real() - lhs[0].imag()*rhs[0].imag() + 
                            lhs[1].real()*rhs[1].real() - lhs[1].imag()*rhs[1].imag() + 
                            lhs[2].real()*rhs[2].real() - lhs[2].imag()*rhs[2].imag() ,
                            lhs[0].real()*rhs[0].imag() + lhs[0].imag()*rhs[0].real() + 
                            lhs[1].real()*rhs[1].imag() + lhs[1].imag()*rhs[1].real() + 
                            lhs[2].real()*rhs[2].imag() + lhs[2].imag()*rhs[2].real() );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Definition of the struct OTCross.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2> struct OTCross {};

//////////////////////////////////////////////////////////////////////
//
// Specializations for TinyVector cross TinyVector
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, unsigned D>
struct OTCross< TinyVector<T1,D> , TinyVector<T2,D> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,D>
  apply(const TinyVector<T1,D>& a, const TinyVector<T2,D>& b)
  {
    TinyVector<Type_t,D> bogusCross(-99999);
    return bogusCross;
  }
};

template<class T1, class T2>
struct OTCross< TinyVector<T1,3> , TinyVector<T2,3> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,3>
  apply(const TinyVector<T1,3>& a, const TinyVector<T2,3>& b)
  {
    TinyVector<Type_t,3> cross;
    cross[0] = a[1]*b[2] - a[2]*b[1];
    cross[1] = a[2]*b[0] - a[0]*b[2];
    cross[2] = a[0]*b[1] - a[1]*b[0];
    return cross;
  }
};

}

/* This has been substituted by OHMMS_META_BINARY_OPERATORS
#define OHMMS_TINYVECTOR_BINARY_OPERATORS(FUNC,TAG)                           \
                                                                              \
template <class T1, class T2, unsigned D>                                     \
struct BinaryReturn< TinyVector<T1,D>, TinyVector<T2,D>, TAG >                \
{                                                                             \
  typedef TinyVector<typename BinaryReturn<T1,T2,TAG>::Type_t, D> Type_t;     \
};                                                                            \
                                                                              \
template <class T1, class T2, unsigned D>                                     \
inline                                                                        \
typename BinaryReturn< TinyVector<T1,D>, TinyVector<T2,D>, TAG >::Type_t      \
FUNC( const TinyVector<T1,D>& v1, const TinyVector<T2,D>& v2 )                \
{                                                                             \
  return OTBinary<TinyVector<T1,D>,TinyVector<T2,D>,TAG>::apply(v1,v2,TAG()); \
}                                                                             \
                                                                              \
template <class T1, class T2, unsigned D>                                     \
struct BinaryReturn< TinyVector<T1,D>, T2, TAG >                              \
{                                                                             \
  typedef TinyVector< typename BinaryReturn<T1,T2,TAG>::Type_t, D > Type_t;   \
};                                                                            \
                                                                              \
template <class T1, class T2, unsigned D>                                     \
struct BinaryReturn< T1, TinyVector<T2,D>, TAG >                              \
{                                                                             \
  typedef TinyVector< typename BinaryReturn<T1,T2,TAG>::Type_t, D > Type_t;   \
};                                                                            \
                                                                              \
template <class T1, class T2, unsigned D>                                     \
inline                                                                        \
typename BinaryReturn< TinyVector<T1,D>, T2, TAG >::Type_t                    \
FUNC( const TinyVector<T1,D>& v1, const T2& x )                               \
{                                                                             \
  return OTBinary<TinyVector<T1,D>,T2,TAG>::apply(v1,x,TAG());                \
}                                                                             \
                                                                              \
template <class T1, class T2, unsigned D>                                     \
inline                                                                        \
typename BinaryReturn< T1, TinyVector<T2,D>, TAG >::Type_t                    \
FUNC( const T1& x, const TinyVector<T2,D>& v2)                                \
{                                                                             \
  return OTBinary<T1,TinyVector<T2,D>,TAG>::apply(x,v2,TAG());                \
}                                                                             \

#define OHMMS_TINYVECTOR_ACCUM_OPERATORS(FUNC,TAG)                            \
                                                                              \
template <class T1, class T2, unsigned D>                                     \
inline TinyVector<T1,D>&                                                      \
FUNC( TinyVector<T1,D>& v1, const TinyVector<T2,D>& v2 )                      \
{                                                                             \
  OTAssign<TinyVector<T1,D>,TinyVector<T2,D>,TAG>::apply(v1,v2,TAG());        \
  return v1;                                                                  \
}                                                                             \
                                                                              \
template <class T1, class T2, unsigned D>                                     \
inline TinyVector<T1,D>&                                                      \
FUNC( TinyVector<T1,D>& v1, const T2& v2 )                                    \
{                                                                             \
  OTAssign<TinyVector<T1,D>,T2,TAG>::apply(v1,v2,TAG());                      \
  return v1;                                                                  \
}                                                                             \

*/
#endif // OHMMS_TINYVECTOR_OPERATORS_H

