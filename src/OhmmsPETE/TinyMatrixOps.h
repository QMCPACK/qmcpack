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


#ifndef OHMMS_TINYMATRIX_OPERATORS_H
#define OHMMS_TINYMATRIX_OPERATORS_H

namespace qmcplusplus
{

///////////////////////////////////////////////////////////////////////
// specialized assignment and element-by-element binary functors for tinymatrix
//
// Generic Assignment operators (defined in OhmmsTinyMeta.h)
// template<class T1, class T2, class OP> struct OTAssign {};
//
// Generic Binary operators (defined in OhmmsTinyMeta.h)
// template<class T1, class T2, class OP> struct OTBinary {};
//
///////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Specializations for TinyMatrixs of arbitrary size.
//////////////////////////////////////////////////////////////////////
template<class T1, class T2, class OP, unsigned D1, unsigned D2>
struct OTAssign< TinyMatrix<T1,D1,D2> , TinyMatrix<T2,D1,D2> , OP >
{
  inline static void
  apply( TinyMatrix<T1,D1,D2>& lhs, const TinyMatrix<T2,D1,D2>& rhs, OP op)
  {
    for (unsigned d=0; d<D1*D2; ++d)
      op(lhs[d] , rhs[d]);
  }
};

template<class T1, class T2, class OP, unsigned D1, unsigned D2>
struct OTAssign< TinyMatrix<T1,D1,D2> , T2 , OP >
{
  inline static void
  apply( TinyMatrix<T1,D1,D2>& lhs, const T2& rhs, OP op )
  {
    for (unsigned d=0; d<D1*D2; ++d)
      op(lhs[d] , rhs);
  }
};


//////////////////////////////////////////////////////////////////////
// Specializations for TinyMatrixs with D=1 x1
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct OTAssign< TinyMatrix<T1,1,1> , TinyMatrix<T2,1,1> , OP >
{
  inline static void
  apply( TinyMatrix<T1,1,1>& lhs, const TinyMatrix<T2,1,1>& rhs, OP op)
  {
    op(lhs[0] , rhs[0] );
  }
};

template<class T1, class T2, class OP>
struct OTAssign< TinyMatrix<T1,1,1> , T2 , OP >
{
  inline static void
  apply( TinyMatrix<T1,1,1>& lhs, const T2& rhs, OP op )
  {
    op(lhs[0] , rhs );
  }
};

//////////////////////////////////////////////////////////////////////
// Specializations for TinyMatrixs with D=1x2 and 2x1
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct OTAssign< TinyMatrix<T1,1,2> , TinyMatrix<T2,1,2> , OP >
{
  inline static void
  apply( TinyMatrix<T1,1,2>& lhs, const TinyMatrix<T2,1,2>& rhs, OP op)
  {
    op(lhs[0] , rhs[0] );
    op(lhs[1] , rhs[1] );
  }
};

template<class T1, class T2, class OP>
struct OTAssign< TinyMatrix<T1,1,2> , T2 , OP >
{
  inline static void
  apply( TinyMatrix<T1,1,2>& lhs, const T2& rhs, OP op )
  {
    op(lhs[0] , rhs);
    op(lhs[1] , rhs);
  }
};

template<class T1, class T2, class OP>
struct OTAssign< TinyMatrix<T1,2,1> , TinyMatrix<T2,2,1> , OP >
{
  inline static void
  apply( TinyMatrix<T1,2,1>& lhs, const TinyMatrix<T2,2,1>& rhs, OP op)
  {
    op(lhs[0] , rhs[0] );
    op(lhs[1] , rhs[1] );
  }
};

template<class T1, class T2, class OP>
struct OTAssign< TinyMatrix<T1,2,1> , T2 , OP >
{
  inline static void
  apply( TinyMatrix<T1,2,1>& lhs, const T2& rhs, OP op )
  {
    op(lhs[0] , rhs);
    op(lhs[1] , rhs);
  }
};
//////////////////////////////////////////////////////////////////////
// Specializations for TinyMatrixs with D=3x1 1x3
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct OTAssign< TinyMatrix<T1,1,3> , TinyMatrix<T2,1,3> , OP >
{
  inline static void
  apply( TinyMatrix<T1,1,3>& lhs, const TinyMatrix<T2,1,3>& rhs, OP op)
  {
    op(lhs[0] , rhs[0] );
    op(lhs[1] , rhs[1] );
    op(lhs[2] , rhs[2] );
  }
};

template<class T1, class T2, class OP>
struct OTAssign< TinyMatrix<T1,1,3> , T2 , OP >
{
  inline static void
  apply( TinyMatrix<T1,1,3>& lhs, const T2& rhs, OP op )
  {
    op(lhs[0] , rhs );
    op(lhs[1] , rhs );
    op(lhs[2] , rhs );
  }
};

template<class T1, class T2, class OP>
struct OTAssign< TinyMatrix<T1,3,1> , TinyMatrix<T2,3,1> , OP >
{
  inline static void
  apply( TinyMatrix<T1,3,1>& lhs, const TinyMatrix<T2,3,1>& rhs, OP op)
  {
    op(lhs[0] , rhs[0] );
    op(lhs[1] , rhs[1] );
    op(lhs[2] , rhs[2] );
  }
};

template<class T1, class T2, class OP>
struct OTAssign< TinyMatrix<T1,3,1> , T2 , OP >
{
  inline static void
  apply( TinyMatrix<T1,3,1>& lhs, const T2& rhs, OP op )
  {
    op(lhs[0] , rhs );
    op(lhs[1] , rhs );
    op(lhs[2] , rhs );
  }
};
//////////////////////////////////////////////////////////////////////
// Specializations for TinyMatrix with D=2x2
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct OTAssign< TinyMatrix<T1,2,2> , TinyMatrix<T2,2,2> , OP >
{
  inline static void
  apply( TinyMatrix<T1,2,2>& lhs, const TinyMatrix<T2,2,2>& rhs, OP op)
  {
    op(lhs[0] , rhs[0] );
    op(lhs[1] , rhs[1] );
    op(lhs[2] , rhs[2] );
    op(lhs[3] , rhs[3] );
  }
};

template<class T1, class T2, class OP>
struct OTAssign< TinyMatrix<T1,2,2> , T2 , OP >
{
  inline static void
  apply( TinyMatrix<T1,2,2>& lhs, const T2& rhs, OP op )
  {
    op(lhs[0] , rhs );
    op(lhs[1] , rhs );
    op(lhs[2] , rhs );
    op(lhs[3] , rhs );
  }
};

//////////////////////////////////////////////////////////////////////
// Specializations for TinyMatrix with D=3x3
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct OTAssign< TinyMatrix<T1,3,3> , TinyMatrix<T2,3,3> , OP >
{
  inline static void
  apply( TinyMatrix<T1,3,3>& lhs, const TinyMatrix<T2,3,3>& rhs, OP op)
  {
    op(lhs[0] , rhs[0] );
    op(lhs[1] , rhs[1] );
    op(lhs[2] , rhs[2] );
    op(lhs[3] , rhs[3] );
    op(lhs[4] , rhs[4] );
    op(lhs[5] , rhs[5] );
    op(lhs[6] , rhs[6] );
    op(lhs[7] , rhs[7] );
    op(lhs[8] , rhs[8] );
  }
};

template<class T1, class T2, class OP>
struct OTAssign< TinyMatrix<T1,3,3> , T2 , OP >
{
  inline static void
  apply( TinyMatrix<T1,3,3>& lhs, const T2& rhs, OP op )
  {
    op(lhs[0] , rhs );
    op(lhs[1] , rhs );
    op(lhs[2] , rhs );
    op(lhs[3] , rhs );
    op(lhs[4] , rhs );
    op(lhs[5] , rhs );
    op(lhs[6] , rhs );
    op(lhs[7] , rhs );
    op(lhs[8] , rhs );
  }
};

//////////////////////////////////////////////////////////////////////
// Specializations for TinyMatrix with D=4x4
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct OTAssign< TinyMatrix<T1,4,4> , TinyMatrix<T2,4,4> , OP >
{
  inline static void
  apply( TinyMatrix<T1,4,4>& lhs, const TinyMatrix<T2,4,4>& rhs, OP op)
  {
    op(lhs[ 0] , rhs[ 0] );
    op(lhs[ 1] , rhs[ 1] );
    op(lhs[ 2] , rhs[ 2] );
    op(lhs[ 3] , rhs[ 3] );
    op(lhs[ 4] , rhs[ 4] );
    op(lhs[ 5] , rhs[ 5] );
    op(lhs[ 6] , rhs[ 6] );
    op(lhs[ 7] , rhs[ 7] );
    op(lhs[ 8] , rhs[ 8] );
    op(lhs[ 9] , rhs[ 9] );
    op(lhs[10] , rhs[10] );
    op(lhs[11] , rhs[11] );
    op(lhs[12] , rhs[12] );
    op(lhs[13] , rhs[13] );
    op(lhs[14] , rhs[14] );
    op(lhs[15] , rhs[15] );
  }
};

template<class T1, class T2, class OP>
struct OTAssign< TinyMatrix<T1,4,4> , T2 , OP >
{
  inline static void
  apply( TinyMatrix<T1,4,4>& lhs, const T2& rhs, OP op )
  {
    op(lhs[ 0] , rhs );
    op(lhs[ 1] , rhs );
    op(lhs[ 2] , rhs );
    op(lhs[ 3] , rhs );
    op(lhs[ 4] , rhs );
    op(lhs[ 5] , rhs );
    op(lhs[ 6] , rhs );
    op(lhs[ 7] , rhs );
    op(lhs[ 8] , rhs );
    op(lhs[ 9] , rhs );
    op(lhs[10] , rhs );
    op(lhs[11] , rhs );
    op(lhs[12] , rhs );
    op(lhs[13] , rhs );
    op(lhs[14] , rhs );
    op(lhs[15] , rhs );
  }
};

///////////////////////////////////////////////////////////////////////
//
// Binary operators
// template<class T1, class T2, class OP> struct OTBinary {};
//
///////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Specializations for TinyMatrixs of arbitrary size.
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP, unsigned D1, unsigned D2>
struct OTBinary< TinyMatrix<T1,D1,D2> , TinyMatrix<T2,D1,D2> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyMatrix<Type_t,D1,D2>
  apply(const TinyMatrix<T1,D1,D2>& lhs, const TinyMatrix<T2,D1,D2>& rhs, OP op)
  {
    TinyMatrix<Type_t,D1,D2> ret;
    for (unsigned d=0; d<D1*D2; ++d)
      ret[d] = op(lhs[d] , rhs[d]);
    return ret;
  }
};

template<class T1, class T2, class OP, unsigned D1, unsigned D2>
struct OTBinary< TinyMatrix<T1,D1,D2> , T2 , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyMatrix<Type_t,D1,D2>
  apply(const TinyMatrix<T1,D1,D2>& lhs, const T2& rhs, OP op)
  {
    TinyMatrix<Type_t,D1,D2> ret;
    for (unsigned d=0 ; d<D1*D2; ++d)
      ret[d] = op(lhs[d] , rhs );
    return ret;
  }
};

template<class T1, class T2, class OP, unsigned D1, unsigned D2>
struct OTBinary< T1, TinyMatrix<T2,D1,D2> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyMatrix<Type_t,D1,D2>
  apply(const T1& lhs, const TinyMatrix<T2,D1,D2>& rhs, OP op)
  {
    TinyMatrix<Type_t,D1,D2> ret;
    for (unsigned d=0; d<D1*D2; ++d)
      ret[d] = op(lhs , rhs[d]);
    return ret;
  }
};

//////////////////////////////////////////////////////////////////////
// Specializations of OTBinary for TinyMatrix with D=1 x 1
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct OTBinary< TinyMatrix<T1,1,1> , TinyMatrix<T2,1,1> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyMatrix<Type_t,1,1>
  apply(const TinyMatrix<T1,1,1>& lhs, const TinyMatrix<T2,1,1>& rhs, OP op)
  {
    return TinyMatrix<Type_t,1,1>( op(lhs[0], rhs[0] ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< TinyMatrix<T1,1,1> , T2 , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyMatrix<Type_t,1,1>
  apply(const TinyMatrix<T1,1,1>& lhs, const T2& rhs, OP op)
  {
    return TinyMatrix<Type_t,1,1>( op(lhs[0], rhs ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< T1, TinyMatrix<T2,1,1> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyMatrix<Type_t,1,1>
  apply(const T1& lhs, const TinyMatrix<T2,1,1>& rhs, OP op)
  {
    return TinyMatrix<Type_t,1,1>( op(lhs, rhs[0] ) );
  }
};

//////////////////////////////////////////////////////////////////////
// Specializations of OTBinary for TinyMatrixs with 1 x 2
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct OTBinary< TinyMatrix<T1,1,2> , TinyMatrix<T2,1,2> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyMatrix<Type_t,1,2>
  apply(const TinyMatrix<T1,1,2>& lhs, const TinyMatrix<T2,1,2>& rhs, OP op)
  {
    return TinyMatrix<Type_t,1,2>( op(lhs[0], rhs[0]),op(lhs[1], rhs[1] ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< TinyMatrix<T1,1,2> , T2 , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyMatrix<Type_t,1,2>
  apply(const TinyMatrix<T1,1,2>& lhs, const T2& rhs, OP op)
  {
    return TinyMatrix<Type_t,1,2>( op(lhs[0], rhs ) ,op(lhs[1], rhs));
  }
};

template<class T1, class T2, class OP>
struct OTBinary< T1, TinyMatrix<T2,1,2> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyMatrix<Type_t,1,2>
  apply(const T1& lhs, const TinyMatrix<T2,1,2>& rhs, OP op)
  {
    return TinyMatrix<Type_t,1,2>( op(lhs, rhs[0] ) ,op(lhs, rhs[1] ) );
  }
};

//////////////////////////////////////////////////////////////////////
// Specializations for TinyMatrixs with D=3x1 1x3
//////////////////////////////////////////////////////////////////////
template<class T1, class T2, class OP>
struct OTBinary< TinyMatrix<T1,1,3> , TinyMatrix<T2,1,3> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyMatrix<Type_t,1,3>
  apply(const TinyMatrix<T1,1,3>& lhs, const TinyMatrix<T2,1,3>& rhs, OP op)
  {
    return TinyMatrix<Type_t,1,3>( op(lhs[0], rhs[0]),
                                   op(lhs[1], rhs[1]),
                                   op(lhs[2], rhs[2]) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< TinyMatrix<T1,1,3> , T2 , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyMatrix<Type_t,1,3>
  apply(const TinyMatrix<T1,1,3>& lhs, const T2& rhs, OP op)
  {
    return TinyMatrix<Type_t,1,3>( op(lhs[0], rhs),
                                   op(lhs[1], rhs),
                                   op(lhs[2], rhs));
  }
};

template<class T1, class T2, class OP>
struct OTBinary< T1, TinyMatrix<T2,1,3> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyMatrix<Type_t,1,3>
  apply(const T1& lhs, const TinyMatrix<T2,1,3>& rhs, OP op)
  {
    return TinyMatrix<Type_t,1,3>( op(lhs, rhs[0]),
                                   op(lhs, rhs[1]),
                                   op(lhs, rhs[2]));
  }
};

//////////////////////////////////////////////////////////////////////
// Specializations of OTBinary for TinyMatrixs with D=2x2
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct OTBinary< TinyMatrix<T1,2,2> , TinyMatrix<T2,2,2> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyMatrix<Type_t,2,2>
  apply(const TinyMatrix<T1,2,2>& lhs, const TinyMatrix<T2,2,2>& rhs, OP op)
  {
    return TinyMatrix<Type_t,2,2>( op(lhs[0], rhs[0] ) ,
                                   op(lhs[1], rhs[1] ) ,
                                   op(lhs[2], rhs[2] ) ,
                                   op(lhs[3], rhs[3] ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< TinyMatrix<T1,2,2> , T2 , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyMatrix<Type_t,2,2>
  apply(const TinyMatrix<T1,2,2>& lhs, const T2& rhs, OP op)
  {
    return TinyMatrix<Type_t,2,2>( op(lhs[0], rhs ) ,
                                   op(lhs[1], rhs ) ,
                                   op(lhs[2], rhs ) ,
                                   op(lhs[3], rhs ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< T1, TinyMatrix<T2,2,2> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static TinyMatrix<Type_t,2,2>
  apply(const T1& lhs, const TinyMatrix<T2,2,2>& rhs, OP op)
  {
    return TinyMatrix<Type_t,2,2>( op( lhs, rhs[0] ) ,
                                   op( lhs, rhs[1] ) ,
                                   op( lhs, rhs[2] ) ,
                                   op( lhs, rhs[3] ) );
  }
};

///////////////////////////////////////////////////////////////
// Specialization for 3x3
///////////////////////////////////////////////////////////////
template<class T1, class T2, class OP>
struct OTBinary< TinyMatrix<T1,3,3> , TinyMatrix<T2,3,3>, OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  typedef TinyMatrix<Type_t, 3, 3> Return_t;
  inline static Return_t
  apply(const TinyMatrix<T1,3,3>& lhs, const TinyMatrix<T2,3,3>& rhs, OP op)
  {
    return Return_t(op(lhs[0],rhs[0]), op(lhs[1],rhs[1]), op(lhs[2],rhs[2]),
                    op(lhs[3],rhs[3]), op(lhs[4],rhs[4]), op(lhs[5],rhs[5]),
                    op(lhs[6],rhs[6]), op(lhs[7],rhs[7]), op(lhs[8],rhs[8]));
  }
};

template<class T1, class T2, class OP>
struct OTBinary< TinyMatrix<T1,3,3> , T2, OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  typedef TinyMatrix<Type_t, 3, 3> Return_t;
  inline static Return_t
  apply(const TinyMatrix<T1,3,3>& lhs, T2 rhs, OP op)
  {
    return Return_t(op(lhs[0],rhs), op(lhs[1],rhs), op(lhs[2],rhs),
                    op(lhs[3],rhs), op(lhs[4],rhs), op(lhs[5],rhs),
                    op(lhs[6],rhs), op(lhs[7],rhs), op(lhs[8],rhs));
  }
};

template<class T1, class T2, class OP>
struct OTBinary< T1, TinyMatrix<T2,3,3>, OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  typedef TinyMatrix<Type_t, 3, 3> Return_t;
  inline static Return_t
  apply(T1 lhs, const TinyMatrix<T2,3,3>& rhs, OP op)
  {
    return Return_t(op(lhs,rhs[0]), op(lhs,rhs[1]), op(lhs,rhs[2]),
                    op(lhs,rhs[3]), op(lhs,rhs[4]), op(lhs,rhs[5]),
                    op(lhs,rhs[6]), op(lhs,rhs[7]), op(lhs,rhs[8]));
  }
};

///////////////////////////////////////////////////////////////
// Specialization for 4x4
///////////////////////////////////////////////////////////////
template<class T1, class T2, class OP>
struct OTBinary< TinyMatrix<T1,4,4> , TinyMatrix<T2,4,4>, OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  typedef TinyMatrix<Type_t, 4, 4> Return_t;
  inline static Return_t
  apply(const TinyMatrix<T1,4,4>& lhs, const TinyMatrix<T2,4,4>& rhs, OP op)
  {
    return Return_t(op(lhs[ 0],rhs[ 0]), op(lhs[ 1],rhs[ 1]), op(lhs[ 2],rhs[ 2]),
                    op(lhs[ 3],rhs[ 3]), op(lhs[ 4],rhs[ 4]), op(lhs[ 5],rhs[ 5]),
                    op(lhs[ 6],rhs[ 6]), op(lhs[ 7],rhs[ 7]), op(lhs[ 8],rhs[ 8]),
                    op(lhs[ 9],rhs[ 9]), op(lhs[10],rhs[10]), op(lhs[11],rhs[11]),
                    op(lhs[12],rhs[12]), op(lhs[13],rhs[13]), op(lhs[14],rhs[14]),
                    op(lhs[15],rhs[15]));
  }
};

template<class T1, class T2, class OP>
struct OTBinary< TinyMatrix<T1,4,4> , T2, OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  typedef TinyMatrix<Type_t, 4, 4> Return_t;
  inline static Return_t
  apply(const TinyMatrix<T1,4,4>& lhs, T2 rhs, OP op)
  {
    return Return_t(op(lhs[ 0],rhs), op(lhs[ 1],rhs), op(lhs[ 2],rhs),
                    op(lhs[ 3],rhs), op(lhs[ 4],rhs), op(lhs[ 5],rhs),
                    op(lhs[ 6],rhs), op(lhs[ 7],rhs), op(lhs[ 8],rhs),
                    op(lhs[ 9],rhs), op(lhs[10],rhs), op(lhs[11],rhs),
                    op(lhs[12],rhs), op(lhs[13],rhs), op(lhs[14],rhs),
                    op(lhs[15],rhs));
  }
};

template<class T1, class T2, class OP>
struct OTBinary< T1, TinyMatrix<T2,4,4>, OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  typedef TinyMatrix<Type_t, 4, 4> Return_t;
  inline static Return_t
  apply(T1 lhs, const TinyMatrix<T2,4,4>& rhs, OP op)
  {
    return Return_t(op(lhs,rhs[ 0]), op(lhs,rhs[ 1]), op(lhs,rhs[ 2]),
                    op(lhs,rhs[ 3]), op(lhs,rhs[ 4]), op(lhs,rhs[ 5]),
                    op(lhs,rhs[ 6]), op(lhs,rhs[ 7]), op(lhs,rhs[ 8]),
                    op(lhs,rhs[ 9]), op(lhs,rhs[10]), op(lhs,rhs[11]),
                    op(lhs,rhs[12]), op(lhs,rhs[13]), op(lhs,rhs[14]),
                    op(lhs,rhs[15]));
  }
};

//////////////////////////////////////////////////////////////////////
//
// ! = Transpose of matrix(D1,D2) to matrix(D2,D1)
//
//////////////////////////////////////////////////////////////////////
template<class T1, unsigned D1, unsigned D2>
inline TinyMatrix<T1,D2,D1>
operator!(const TinyMatrix<T1,D1,D2>& a )
{
  TinyMatrix<T1,D1,D2> res;
  for(int i=0; i<D1; i++)
    for(int j=0; j<D2; j++)
      res(i,j) = a(j,i);
  return res;
}

//////////////////////////////////////////////////////////////////////
// ! specialization for 3x3
//////////////////////////////////////////////////////////////////////
template<class T1>
inline TinyMatrix<T1,3,3>
operator!(const TinyMatrix<T1,3,3>& a )
{
  return TinyMatrix<T1,3,3>(a(0,0),a(1,0),a(2,0),
                            a(0,1),a(1,1),a(2,1),
                            a(0,2),a(1,2),a(2,2));
}

//////////////////////////////////////////////////////////////////////
// ! specialization for 4x4
//////////////////////////////////////////////////////////////////////
template<class T1>
inline TinyMatrix<T1,4,4>
operator!(const TinyMatrix<T1,4,4>& a )
{
  return TinyMatrix<T1,4,4>(a(0,0),a(1,0),a(2,0),a(3,0),
                            a(0,1),a(1,1),a(2,1),a(3,1),
                            a(0,2),a(1,2),a(2,2),a(3,2),
                            a(0,3),a(1,3),a(2,3),a(3,3));
}

//////////////////////////////////////////////////////////////////////
//
// template<class T1, class T2> struct OTDot {};
// Specializations for TinyMatrix x TinyMatrix matrix multiplication
// Matrix(D1,D2)* Matrix(D2,D3)  = Matrix(D1,D3)
//
//////////////////////////////////////////////////////////////////////
template<class T1, class T2, unsigned D1, unsigned D2, unsigned D3>
struct OTDot< TinyMatrix<T1,D1,D2> , TinyMatrix<T2,D2,D3> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  typedef TinyMatrix<Type_t, D1, D3> Return_t;
  inline static Return_t
  apply(const TinyMatrix<T1,D1,D2>& lhs, const TinyMatrix<T2,D2,D3>& rhs)
  {
    Return_t res;
    for(int i=0; i<D1; i++)
      for(int j=0; j<D3; j++)
      {
        Type_t tmp = 0.0e0;
        for(int k=0; k<D2; k++)
          tmp += lhs(i,k)*rhs(k,j);
        res(i,j) = tmp;
      }
    return res;
  }
};

template<class T1, class T2>
struct OTDot< TinyMatrix<T1,3,3> , TinyMatrix<T2,3,3> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  typedef TinyMatrix<Type_t, 3, 3> Return_t;
  inline static Return_t
  apply(const TinyMatrix<T1,3,3>& lhs, const TinyMatrix<T2,3,3>& rhs)
  {
    return Return_t(lhs(0,0)*rhs(0,0)+lhs(0,1)*rhs(1,0)+lhs(0,2)*rhs(2,0),
                    lhs(0,0)*rhs(0,1)+lhs(0,1)*rhs(1,1)+lhs(0,2)*rhs(2,1),
                    lhs(0,0)*rhs(0,2)+lhs(0,1)*rhs(1,2)+lhs(0,2)*rhs(2,2),
                    lhs(1,0)*rhs(0,0)+lhs(1,1)*rhs(1,0)+lhs(1,2)*rhs(2,0),
                    lhs(1,0)*rhs(0,1)+lhs(1,1)*rhs(1,1)+lhs(1,2)*rhs(2,1),
                    lhs(1,0)*rhs(0,2)+lhs(1,1)*rhs(1,2)+lhs(1,2)*rhs(2,2),
                    lhs(2,0)*rhs(0,0)+lhs(2,1)*rhs(1,0)+lhs(2,2)*rhs(2,0),
                    lhs(2,0)*rhs(0,1)+lhs(2,1)*rhs(1,1)+lhs(2,2)*rhs(2,1),
                    lhs(2,0)*rhs(0,2)+lhs(2,1)*rhs(1,2)+lhs(2,2)*rhs(2,2));
  }
};

//////////////////////////////////////////////////////////////////////
//
// template<class T1, class T2> struct OTDot {};
// Specializations for TinyMatrix x TinyVector multiplication
//
//////////////////////////////////////////////////////////////////////
template<class T1, class T2, unsigned D1, unsigned D2>
struct OTDot< TinyMatrix<T1,D1,D2> , TinyVector<T2,D2> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  typedef TinyVector<Type_t, D1> Return_t;
  inline static Return_t
  apply(const TinyMatrix<T1,D1,D2>& lhs, const TinyVector<T2,D2>& rhs)
  {
    Return_t res;
    for(int i=0; i<D1; i++)
    {
      res(i) = lhs(i,0)*rhs(0);
      for(int j=1; j<D2; j++)
        res(i)+= lhs(i,j)*rhs(j);
    }
    return res;
  }
};

template<class T1, class T2, unsigned D1, unsigned D2>
struct OTDot< TinyVector<T1,D1>, TinyMatrix<T2,D1,D2> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  typedef TinyVector<Type_t, D2> Return_t;
  inline static Return_t
  apply(const TinyVector<T1,D1>& lhs, const TinyMatrix<T2,D1,D2>& rhs)
  {
    Return_t res;
    for(int i=0; i<D2; i++)
    {
      res(i) = lhs(0)*rhs(0,i);
      for(int j=1; j<D1; j++)
        res(i)+= lhs(j)*rhs(j,i);
    }
    return res;
  }
};


//////////////////////////////////////////////////////////////////////
// Specializations for TinyMatrixs with D=3x3
//////////////////////////////////////////////////////////////////////
template<class T1, class T2>
struct OTDot< TinyMatrix<T1,3,3> , TinyVector<T2,3> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  typedef TinyVector<Type_t, 3> Return_t;
  inline static Return_t
  apply(const TinyMatrix<T1,3,3>& lhs, const TinyVector<T2,3>& rhs)
  {
    return Return_t(lhs(0,0)*rhs(0) + lhs(0,1)*rhs(1) + lhs(0,2)*rhs(2),
                    lhs(1,0)*rhs(0) + lhs(1,1)*rhs(1) + lhs(1,2)*rhs(2),
                    lhs(2,0)*rhs(0) + lhs(2,1)*rhs(1) + lhs(2,2)*rhs(2));
  }
};

template<class T1, class T2>
struct OTDot< TinyVector<T1,3>, TinyMatrix<T2,3,3> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  typedef TinyVector<Type_t, 3> Return_t;
  inline static Return_t
  apply(const TinyVector<T1,3>& lhs, const TinyMatrix<T2,3,3>& rhs)
  {
    return Return_t(lhs(0)*rhs(0,0) + lhs(1)*rhs(1,0) + lhs(2)*rhs(2,0),
                    lhs(0)*rhs(0,1) + lhs(1)*rhs(1,1) + lhs(2)*rhs(2,1),
                    lhs(0)*rhs(0,2) + lhs(1)*rhs(1,2) + lhs(2)*rhs(2,2));
  }
};

//////////////////////////////////////////////////////////////////////
// Specializations for TinyMatrixs with D=4x4
//////////////////////////////////////////////////////////////////////
template<class T1, class T2>
struct OTDot< TinyMatrix<T1,4,4> , TinyVector<T2,4> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  typedef TinyVector<Type_t, 4> Return_t;
  inline static Return_t
  apply(const TinyMatrix<T1,4,4>& lhs, const TinyVector<T2,4>& rhs)
  {
    return Return_t(lhs(0,0)*rhs(0)+lhs(0,1)*rhs(1)+lhs(0,2)*rhs(2)+lhs(0,3)*rhs(3),
                    lhs(1,0)*rhs(0)+lhs(1,1)*rhs(1)+lhs(1,2)*rhs(2)+lhs(1,3)*rhs(3),
                    lhs(2,0)*rhs(0)+lhs(2,1)*rhs(1)+lhs(2,2)*rhs(2)+lhs(2,3)*rhs(3),
                    lhs(3,0)*rhs(0)+lhs(3,1)*rhs(1)+lhs(3,2)*rhs(2)+lhs(3,3)*rhs(3));
  }
};

template<class T1, class T2>
struct OTDot< TinyVector<T1,4>, TinyMatrix<T2,4,4> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  typedef TinyVector<Type_t, 4> Return_t;
  inline static Return_t
  apply(const TinyVector<T1,4>& lhs, const TinyMatrix<T2,4,4>& rhs)
  {
    return Return_t(lhs(0)*rhs(0,0) + lhs(1)*rhs(1,0) + lhs(2)*rhs(2,0) + lhs(3)*rhs(3,0),
                    lhs(0)*rhs(0,1) + lhs(1)*rhs(1,1) + lhs(2)*rhs(2,1) + lhs(3)*rhs(3,1),
                    lhs(0)*rhs(0,2) + lhs(1)*rhs(1,2) + lhs(2)*rhs(2,2) + lhs(3)*rhs(3,2),
                    lhs(0)*rhs(0,3) + lhs(1)*rhs(1,3) + lhs(2)*rhs(2,3) + lhs(3)*rhs(3,3));
  }
};


template<class T1, class T2>
struct OTDot< TinyMatrix<T1,4,4> , TinyMatrix<T2,4,4> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  typedef TinyMatrix<Type_t, 4, 4> Return_t;
  inline static Return_t
  apply(const TinyMatrix<T1,4,4>& lhs, const TinyMatrix<T2,4,4>& rhs)
  {
    return Return_t(lhs(0,0)*rhs(0,0)+lhs(0,1)*rhs(1,0)+lhs(0,2)*rhs(2,0)+lhs(0,3)*rhs(3,0),
                    lhs(0,0)*rhs(0,1)+lhs(0,1)*rhs(1,1)+lhs(0,2)*rhs(2,1)+lhs(0,3)*rhs(3,1),
                    lhs(0,0)*rhs(0,2)+lhs(0,1)*rhs(1,2)+lhs(0,2)*rhs(2,2)+lhs(0,3)*rhs(3,2),
                    lhs(0,0)*rhs(0,3)+lhs(0,1)*rhs(1,3)+lhs(0,2)*rhs(2,3)+lhs(0,3)*rhs(3,3),
                    lhs(1,0)*rhs(0,0)+lhs(1,1)*rhs(1,0)+lhs(1,2)*rhs(2,0)+lhs(1,3)*rhs(3,0),
                    lhs(1,0)*rhs(0,1)+lhs(1,1)*rhs(1,1)+lhs(1,2)*rhs(2,1)+lhs(1,3)*rhs(3,1),
                    lhs(1,0)*rhs(0,2)+lhs(1,1)*rhs(1,2)+lhs(1,2)*rhs(2,2)+lhs(1,3)*rhs(3,2),
                    lhs(1,0)*rhs(0,3)+lhs(1,1)*rhs(1,3)+lhs(1,2)*rhs(2,3)+lhs(1,3)*rhs(3,3),
                    lhs(2,0)*rhs(0,0)+lhs(2,1)*rhs(1,0)+lhs(2,2)*rhs(2,0)+lhs(2,3)*rhs(3,0),
                    lhs(2,0)*rhs(0,1)+lhs(2,1)*rhs(1,1)+lhs(2,2)*rhs(2,1)+lhs(2,3)*rhs(3,1),
                    lhs(2,0)*rhs(0,2)+lhs(2,1)*rhs(1,2)+lhs(2,2)*rhs(2,2)+lhs(2,3)*rhs(3,2),
                    lhs(2,0)*rhs(0,3)+lhs(2,1)*rhs(1,3)+lhs(2,2)*rhs(2,3)+lhs(2,3)*rhs(3,3),
                    lhs(3,0)*rhs(0,0)+lhs(3,1)*rhs(1,0)+lhs(3,2)*rhs(2,0)+lhs(3,3)*rhs(3,0),
                    lhs(3,0)*rhs(0,1)+lhs(3,1)*rhs(1,1)+lhs(3,2)*rhs(2,1)+lhs(3,3)*rhs(3,1),
                    lhs(3,0)*rhs(0,2)+lhs(3,1)*rhs(1,2)+lhs(3,2)*rhs(2,2)+lhs(3,3)*rhs(3,2),
                    lhs(3,0)*rhs(0,3)+lhs(3,1)*rhs(1,3)+lhs(3,2)*rhs(2,3)+lhs(3,3)*rhs(3,3));
  }
};

}
#endif // OHMMS_TINYVECTOR_OPERATORS_H

