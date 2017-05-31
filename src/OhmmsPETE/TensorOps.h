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


#ifndef OHMMS_TENSOR_OPERATORS_H
#define OHMMS_TENSOR_OPERATORS_H

/*** Tenor operators.  Generic operators are specialized for 1,2 and 3 D
 */
namespace qmcplusplus
{

template<class T1, class T2, class OP, unsigned D>
struct  OTAssign< Tensor<T1,D> , Tensor<T2,D> , OP >
{
  inline static void
  apply( Tensor<T1,D>& lhs, const Tensor<T2,D>& rhs, OP op)
  {
    for (unsigned d=0; d<D*D; ++d)
      op(lhs[d] , rhs[d]);
  }
};

template<class T1, class T2, class OP, unsigned D>
struct  OTAssign< Tensor<T1,D> , T2 , OP >
{
  inline static void
  apply( Tensor<T1,D>& lhs, T2 rhs, OP op)
  {
    for (unsigned d=0; d<D*D; ++d)
      op(lhs[d] , rhs);
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations for Tensors with D=1.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct  OTAssign< Tensor<T1,1> , Tensor<T2,1> , OP >
{
  inline static void
  apply( Tensor<T1,1>& lhs, const Tensor<T2,1>& rhs, OP op)
  {
    op(lhs[0] , rhs[0] );
  }
};

template<class T1, class T2, class OP>
struct  OTAssign< Tensor<T1,1> , T2 , OP >
{
  inline static void
  apply( Tensor<T1,1>& lhs, T2 rhs, OP op )
  {
    op(lhs[0] , rhs );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations for Tensors with D=2.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct  OTAssign< Tensor<T1,2> , Tensor<T2,2> , OP >
{
  inline static void
  apply( Tensor<T1,2>& lhs, const Tensor<T2,2>& rhs, OP op)
  {
    op(lhs[0] , rhs[0] );
    op(lhs[1] , rhs[1] );
    op(lhs[2] , rhs[2] );
    op(lhs[3] , rhs[3] );
  }
};

template<class T1, class T2, class OP>
struct  OTAssign< Tensor<T1,2> , T2 , OP >
{
  inline static void
  apply( Tensor<T1,2>& lhs, T2 rhs, OP op )
  {
    op(lhs[0] , rhs );
    op(lhs[1] , rhs );
    op(lhs[2] , rhs );
    op(lhs[3] , rhs );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations for Tensors with D=3.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct  OTAssign< Tensor<T1,3> , Tensor<T2,3> , OP >
{
  inline static void
  apply( Tensor<T1,3>& lhs, const Tensor<T2,3>& rhs, OP op)
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
struct  OTAssign< Tensor<T1,3> , T2 , OP >
{
  inline static void
  apply( Tensor<T1,3>& lhs, T2 rhs, OP op )
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
//
// The default definitions for SymTensors of arbitrary size.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP, unsigned D>
struct  OTAssign< SymTensor<T1,D> , SymTensor<T2,D> , OP >
{
  inline static void
  apply( SymTensor<T1,D>& lhs, const SymTensor<T2,D>& rhs, OP op)
  {
    for (unsigned d=0; d<D*(D+1)/2; ++d)
      op(lhs[d] , rhs[d]);
  }
};

template<class T1, class T2, class OP, unsigned D>
struct  OTAssign< SymTensor<T1,D> , T2 , OP >
{
  inline static void
  apply( SymTensor<T1,D>& lhs, T2 rhs, OP op )
  {
    for (unsigned d=0; d<D*(D+1)/2; ++d)
      op(lhs[d] , rhs);
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations for SymTensors with D=1.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct  OTAssign< SymTensor<T1,1> , SymTensor<T2,1> , OP >
{
  inline static void
  apply( SymTensor<T1,1>& lhs, const SymTensor<T2,1>& rhs, OP op)
  {
    op(lhs[0] , rhs[0] );
  }
};

template<class T1, class T2, class OP>
struct  OTAssign< SymTensor<T1,1> , T2 , OP>
{
  inline static void
  apply( SymTensor<T1,1>& lhs, T2 rhs, OP op )
  {
    op(lhs[0] , rhs );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations for SymTensors with D=2.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct  OTAssign< SymTensor<T1,2> , SymTensor<T2,2> , OP >
{
  inline static void
  apply( SymTensor<T1,2>& lhs, const SymTensor<T2,2>& rhs, OP op)
  {
    op(lhs[0] , rhs[0] );
    op(lhs[1] , rhs[1] );
    op(lhs[2] , rhs[2] );
  }
};

template<class T1, class T2, class OP>
struct  OTAssign< SymTensor<T1,2> , T2 , OP >
{
  inline static void
  apply( SymTensor<T1,2>& lhs, T2 rhs, OP op )
  {
    op(lhs[0] , rhs );
    op(lhs[1] , rhs );
    op(lhs[2] , rhs );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations for SymTensors with D=3.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct  OTAssign< SymTensor<T1,3> , SymTensor<T2,3> , OP >
{
  inline static void
  apply( SymTensor<T1,3>& lhs, const SymTensor<T2,3>& rhs, OP op)
  {
    op(lhs[0] , rhs[0] );
    op(lhs[1] , rhs[1] );
    op(lhs[2] , rhs[2] );
    op(lhs[3] , rhs[3] );
    op(lhs[4] , rhs[4] );
    op(lhs[5] , rhs[5] );
  }
};

template<class T1, class T2, class OP>
struct  OTAssign< SymTensor<T1,3> , T2 , OP >
{
  inline static void
  apply( SymTensor<T1,3>& lhs, T2 rhs, OP op )
  {
    op(lhs[0] , rhs );
    op(lhs[1] , rhs );
    op(lhs[2] , rhs );
    op(lhs[3] , rhs );
    op(lhs[4] , rhs );
    op(lhs[5] , rhs );
  }
};

//////////////////////////////////////////////////////////////////////
//
// The default definitions for AntiSymTensors of arbitrary size.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP, unsigned D>
struct  OTAssign< AntiSymTensor<T1,D> , AntiSymTensor<T2,D> , OP >
{
  inline static void
  apply( AntiSymTensor<T1,D>& lhs, const AntiSymTensor<T2,D>& rhs, OP op)
  {
    for (unsigned d=0; d<D*(D-1)/2; ++d)
      op(lhs[d] , rhs[d]);
  }
};

template<class T1, class T2, class OP, unsigned D>
struct  OTAssign< AntiSymTensor<T1,D> , T2 , OP >
{
  inline static void
  apply( AntiSymTensor<T1,D>& lhs, T2 rhs, OP op )
  {
    for (unsigned d=0; d<D*(D-1)/2; ++d)
      op(lhs[d] , rhs);
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations for AntiSymTensors with D=1.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct  OTAssign< AntiSymTensor<T1,1> , AntiSymTensor<T2,1> , OP >
{
  inline static void
  apply( AntiSymTensor<T1,1>& lhs, const AntiSymTensor<T2,1>& rhs, OP op)
  {
  }
};

template<class T1, class T2, class OP>
struct  OTAssign< AntiSymTensor<T1,1> , T2 , OP >
{
  inline static void
  apply( AntiSymTensor<T1,1>& lhs, T2 rhs, OP op )
  {
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations for AntiSymTensors with D=2.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct  OTAssign< AntiSymTensor<T1,2> , AntiSymTensor<T2,2> , OP >
{
  inline static void
  apply( AntiSymTensor<T1,2>& lhs, const AntiSymTensor<T2,2>& rhs, OP op)
  {
    op(lhs[0] , rhs[0] );
  }
};

template<class T1, class T2, class OP>
struct  OTAssign< AntiSymTensor<T1,2> , T2 , OP >
{
  inline static void
  apply( AntiSymTensor<T1,2>& lhs, T2 rhs, OP op )
  {
    op(lhs[0] , rhs );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations for AntiSymTensors with D=3.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct  OTAssign< AntiSymTensor<T1,3> , AntiSymTensor<T2,3> , OP >
{
  inline static void
  apply( AntiSymTensor<T1,3>& lhs, const AntiSymTensor<T2,3>& rhs, OP op)
  {
    op(lhs[0] , rhs[0] );
    op(lhs[1] , rhs[1] );
    op(lhs[2] , rhs[2] );
  }
};

template<class T1, class T2, class OP>
struct  OTAssign< AntiSymTensor<T1,3> , T2 , OP >
{
  inline static void
  apply( AntiSymTensor<T1,3>& lhs, T2 rhs, OP op )
  {
    op(lhs[0] , rhs );
    op(lhs[1] , rhs );
    op(lhs[2] , rhs );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations for Tensors of arbitrary size.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP, unsigned D>
struct OTBinary< Tensor<T1,D> , Tensor<T2,D> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static Tensor<Type_t,D>
  apply(const Tensor<T1,D>& lhs, const Tensor<T2,D>& rhs, OP op)
  {
    Tensor<Type_t,D> ret;
    for (unsigned d=0; d<D*D; ++d)
      ret[d] = op( lhs[d] , rhs[d]);
    return ret;
  }
};

template<class T1, class T2, class OP, unsigned D>
struct OTBinary< Tensor<T1,D> , T2 , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static Tensor<Type_t,D>
  apply(const Tensor<T1,D>& lhs, T2 rhs, OP op)
  {
    Tensor<Type_t,D> ret;
    for (unsigned d=0; d<D*D; ++d)
      ret[d] = op(lhs[d] , rhs );
    return ret;
  }
};

template<class T1, class T2, class OP, unsigned D>
struct OTBinary< T1, Tensor<T2,D> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static Tensor<Type_t,D>
  apply(T1 lhs, const Tensor<T2,D>& rhs, OP op)
  {
    Tensor<Type_t,D> ret;
    for (unsigned d=0; d<D*D; ++d)
      ret[d] = op(lhs , rhs[d]);
    return ret;
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations of OTBinary for Tensors with D=1.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct OTBinary< Tensor<T1,1> , Tensor<T2,1> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static Tensor<Type_t,1>
  apply(const Tensor<T1,1>& lhs, const Tensor<T2,1>& rhs, OP op)
  {
    return Tensor<Type_t,1>( op( lhs[0], rhs[0] ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< Tensor<T1,1> , T2 , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static Tensor<Type_t,1>
  apply(const Tensor<T1,1>& lhs, T2 rhs, OP op)
  {
    return Tensor<Type_t,1>( op(lhs[0], rhs ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< T1, Tensor<T2,1> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static Tensor<Type_t,1>
  apply(T1 lhs, const Tensor<T2,1>& rhs, OP op)
  {
    return Tensor<Type_t,1>( op( lhs, rhs[0] ) );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations of OTBinary for Tensors with D=2.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct OTBinary< Tensor<T1,2> , Tensor<T2,2> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static Tensor<Type_t,2>
  apply(const Tensor<T1,2>& lhs, const Tensor<T2,2>& rhs, OP op)
  {
    return Tensor<Type_t,2>( op( lhs[0], rhs[0] ) ,
                             op( lhs[1], rhs[1] ) ,
                             op( lhs[2], rhs[2] ) ,
                             op( lhs[3], rhs[3] ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< Tensor<T1,2> , T2 , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static Tensor<Type_t,2>
  apply(const Tensor<T1,2>& lhs, T2 rhs, OP op)
  {
    return Tensor<Type_t,2>( op( lhs[0], rhs ) ,
                             op( lhs[1], rhs ) ,
                             op( lhs[2], rhs ) ,
                             op( lhs[3], rhs ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< T1, Tensor<T2,2> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static Tensor<Type_t,2>
  apply(T1 lhs, const Tensor<T2,2>& rhs, OP op)
  {
    return Tensor<Type_t,2>( op( lhs, rhs[0] ) ,
                             op( lhs, rhs[1] ) ,
                             op( lhs, rhs[2] ) ,
                             op( lhs, rhs[3] ) );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations of OTBinary for Tensors with D=3.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct OTBinary< Tensor<T1,3> , Tensor<T2,3> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static Tensor<Type_t,3>
  apply(const Tensor<T1,3>& lhs, const Tensor<T2,3>& rhs, OP op)
  {
    return Tensor<Type_t,3>( op( lhs[0], rhs[0] ) ,
                             op( lhs[1], rhs[1] ) ,
                             op( lhs[2], rhs[2] ) ,
                             op( lhs[3], rhs[3] ) ,
                             op( lhs[4], rhs[4] ) ,
                             op( lhs[5], rhs[5] ) ,
                             op( lhs[6], rhs[6] ) ,
                             op( lhs[7], rhs[7] ) ,
                             op( lhs[8], rhs[8] ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< Tensor<T1,3> , T2 , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static Tensor<Type_t,3>
  apply(const Tensor<T1,3>& lhs, T2 rhs, OP op)
  {
    return Tensor<Type_t,3>( op( lhs[0], rhs ) ,
                             op( lhs[1], rhs ) ,
                             op( lhs[2], rhs ) ,
                             op( lhs[3], rhs ) ,
                             op( lhs[4], rhs ) ,
                             op( lhs[5], rhs ) ,
                             op( lhs[6], rhs ) ,
                             op( lhs[7], rhs ) ,
                             op( lhs[8], rhs ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< T1, Tensor<T2,3> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static Tensor<Type_t,3>
  apply(T1 lhs, const Tensor<T2,3>& rhs, OP op)
  {
    return Tensor<Type_t,3>( op( lhs, rhs[0] ) ,
                             op( lhs, rhs[1] ) ,
                             op( lhs, rhs[2] ) ,
                             op( lhs, rhs[3] ) ,
                             op( lhs, rhs[4] ) ,
                             op( lhs, rhs[5] ) ,
                             op( lhs, rhs[6] ) ,
                             op( lhs, rhs[7] ) ,
                             op( lhs, rhs[8] ) );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations for SymTensors of arbitrary size.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP, unsigned D>
struct OTBinary< SymTensor<T1,D> , SymTensor<T2,D> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static SymTensor<Type_t,D>
  apply(const SymTensor<T1,D>& lhs, const SymTensor<T2,D>& rhs, OP op)
  {
    SymTensor<Type_t,D> ret;
    for (unsigned d=0; d<D*(D+1)/2; ++d)
      ret[d] = op(lhs[d] , rhs[d]);
    return ret;
  }
};

template<class T1, class T2, class OP, unsigned D>
struct OTBinary< SymTensor<T1,D> , T2 , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static SymTensor<Type_t,D>
  apply(const SymTensor<T1,D>& lhs, T2 rhs, OP op)
  {
    SymTensor<Type_t,D> ret;
    for (unsigned d=0; d<D*(D+1)/2; ++d)
      ret[d] = op( lhs[d] , rhs );
    return ret;
  }
};

template<class T1, class T2, class OP, unsigned D>
struct OTBinary< T1, SymTensor<T2,D> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static SymTensor<Type_t,D>
  apply(T1 lhs, const SymTensor<T2,D>& rhs, OP op)
  {
    SymTensor<Type_t,D> ret;
    for (unsigned d=0; d<D*(D+1)/2; ++d)
      ret[d] = op( lhs , rhs[d]);
    return ret;
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations of OTBinary for SymTensors with D=1.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct OTBinary< SymTensor<T1,1> , SymTensor<T2,1> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static SymTensor<Type_t,1>
  apply(const SymTensor<T1,1>& lhs, const SymTensor<T2,1>& rhs, OP op)
  {
    return SymTensor<Type_t,1>( op( lhs[0], rhs[0] ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< SymTensor<T1,1> , T2 , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static SymTensor<Type_t,1>
  apply(const SymTensor<T1,1>& lhs, T2 rhs, OP op)
  {
    return SymTensor<Type_t,1>( op( lhs[0], rhs ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< T1, SymTensor<T2,1> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static SymTensor<Type_t,1>
  apply(T1 lhs, const SymTensor<T2,1>& rhs, OP op)
  {
    return SymTensor<Type_t,1>( op( lhs, rhs[0] ) );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations of OTBinary for SymTensors with D=2.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct OTBinary< SymTensor<T1,2> , SymTensor<T2,2> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static SymTensor<Type_t,2>
  apply(const SymTensor<T1,2>& lhs, const SymTensor<T2,2>& rhs, OP op)
  {
    return SymTensor<Type_t,2>( op( lhs[0], rhs[0] ) ,
                                op( lhs[1], rhs[1] ) ,
                                op( lhs[2], rhs[2] ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< SymTensor<T1,2> , T2 , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static SymTensor<Type_t,2>
  apply(const SymTensor<T1,2>& lhs, T2 rhs, OP op)
  {
    return SymTensor<Type_t,2>( op( lhs[0], rhs ) ,
                                op( lhs[1], rhs ) ,
                                op( lhs[2], rhs ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< T1, SymTensor<T2,2> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static SymTensor<Type_t,2>
  apply(T1 lhs, const SymTensor<T2,2>& rhs, OP op)
  {
    return SymTensor<Type_t,2>( op( lhs, rhs[0] ) ,
                                op( lhs, rhs[1] ) ,
                                op( lhs, rhs[2] ) );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations of OTBinary for SymTensors with D=3.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct OTBinary< SymTensor<T1,3> , SymTensor<T2,3> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static SymTensor<Type_t,3>
  apply(const SymTensor<T1,3>& lhs, const SymTensor<T2,3>& rhs, OP op)
  {
    return SymTensor<Type_t,3>( op( lhs[0], rhs[0] ) ,
                                op( lhs[1], rhs[1] ) ,
                                op( lhs[2], rhs[2] ) ,
                                op( lhs[3], rhs[3] ) ,
                                op( lhs[4], rhs[4] ) ,
                                op( lhs[5], rhs[5] ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< SymTensor<T1,3> , T2 , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static SymTensor<Type_t,3>
  apply(const SymTensor<T1,3>& lhs, T2 rhs, OP op)
  {
    return SymTensor<Type_t,3>( op( lhs[0], rhs ) ,
                                op( lhs[1], rhs ) ,
                                op( lhs[2], rhs ) ,
                                op( lhs[3], rhs ) ,
                                op( lhs[4], rhs ) ,
                                op( lhs[5], rhs ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< T1, SymTensor<T2,3> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static SymTensor<Type_t,3>
  apply(T1 lhs, const SymTensor<T2,3>& rhs, OP op)
  {
    return SymTensor<Type_t,3>( op( lhs, rhs[0] ) ,
                                op( lhs, rhs[1] ) ,
                                op( lhs, rhs[2] ) ,
                                op( lhs, rhs[3] ) ,
                                op( lhs, rhs[4] ) ,
                                op( lhs, rhs[5] ) );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specialization for SymTensor OP Tensor of arbitrary size.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP, unsigned D>
struct OTBinary< SymTensor<T1,D>, Tensor<T2,D>, OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static Tensor<Type_t,D>
  apply(const SymTensor<T1,D> &lhs, const Tensor<T2,D> &rhs, OP op)
  {
    Tensor<Type_t,D> ret;
    for (unsigned i = 0; i < D; i++)
      for (unsigned j = 0; j < D; j++)
        ret(i, j) = op(lhs(i, j), rhs(i, j));
    return ret;
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specialization for Tensor OP SymTensor of arbitrary size.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP, unsigned D>
struct OTBinary< Tensor<T1,D>, SymTensor<T2,D>, OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static Tensor<Type_t,D>
  apply(const Tensor<T1,D> &lhs, const SymTensor<T2,D> &rhs, OP op)
  {
    Tensor<Type_t,D> ret;
    for (unsigned i = 0; i < D; i++)
      for (unsigned j = 0; j < D; j++)
        ret(i, j) = op( lhs(i, j), rhs(i, j));
    return ret;
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations for AntiSymTensors of arbitrary size.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP, unsigned D>
struct OTBinary< AntiSymTensor<T1,D> , AntiSymTensor<T2,D> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static AntiSymTensor<Type_t,D>
  apply(const AntiSymTensor<T1,D>& lhs, const AntiSymTensor<T2,D>& rhs, OP op)
  {
    AntiSymTensor<Type_t,D> ret;
    for (unsigned d=0; d<D*(D-1)/2; ++d)
      ret[d] = op(lhs[d] , rhs[d]);
    return ret;
  }
};

template<class T1, class T2, class OP, unsigned D>
struct OTBinary< AntiSymTensor<T1,D> , T2 , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static AntiSymTensor<Type_t,D>
  apply(const AntiSymTensor<T1,D>& lhs, T2 rhs, OP op)
  {
    AntiSymTensor<Type_t,D> ret;
    for (unsigned d=0; d<D*(D-1)/2; ++d)
      ret[d] = op( lhs[d] , rhs );
    return ret;
  }
};

template<class T1, class T2, class OP, unsigned D>
struct OTBinary< T1, AntiSymTensor<T2,D> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static AntiSymTensor<Type_t,D>
  apply(T1 lhs, const AntiSymTensor<T2,D>& rhs, OP op)
  {
    AntiSymTensor<Type_t,D> ret;
    for (unsigned d=0; d<D*(D-1)/2; ++d)
      ret[d] = op( lhs , rhs[d]);
    return ret;
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations of OTBinary for AntiSymTensors with D=1.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct OTBinary< AntiSymTensor<T1,1> , AntiSymTensor<T2,1> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static AntiSymTensor<Type_t,1>
  apply(const AntiSymTensor<T1,1>& lhs, const AntiSymTensor<T2,1>& rhs, OP op)
  {
    return AntiSymTensor<Type_t,1>( AntiSymTensor<Type_t,1>::DontInitialize() );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< AntiSymTensor<T1,1> , T2 , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static AntiSymTensor<Type_t,1>
  apply(const AntiSymTensor<T1,1>& lhs, T2 rhs, OP op)
  {
    return AntiSymTensor<Type_t,1>( AntiSymTensor<Type_t,1>::DontInitialize() );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< T1, AntiSymTensor<T2,1> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static AntiSymTensor<Type_t,1>
  apply(T1 lhs, const AntiSymTensor<T2,1>& rhs, OP op)
  {
    return AntiSymTensor<Type_t,1>( AntiSymTensor<Type_t,1>::DontInitialize() );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations of OTBinary for AntiSymTensors with D=2.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct OTBinary< AntiSymTensor<T1,2> , AntiSymTensor<T2,2> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static AntiSymTensor<Type_t,2>
  apply(const AntiSymTensor<T1,2>& lhs, const AntiSymTensor<T2,2>& rhs, OP op)
  {
    return AntiSymTensor<Type_t,2>( op( lhs[0], rhs[0] ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< AntiSymTensor<T1,2> , T2 , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static AntiSymTensor<Type_t,2>
  apply(const AntiSymTensor<T1,2>& lhs, T2 rhs, OP op)
  {
    return AntiSymTensor<Type_t,2>( op( lhs[0], rhs ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< T1, AntiSymTensor<T2,2> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static AntiSymTensor<Type_t,2>
  apply(T1 lhs, const AntiSymTensor<T2,2>& rhs, OP op)
  {
    return AntiSymTensor<Type_t,2>( op( lhs, rhs[0] ) );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations of OTBinary for AntiSymTensors with D=3.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, class OP>
struct OTBinary< AntiSymTensor<T1,3> , AntiSymTensor<T2,3> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static AntiSymTensor<Type_t,3>
  apply(const AntiSymTensor<T1,3>& lhs, const AntiSymTensor<T2,3>& rhs, OP op)
  {
    return AntiSymTensor<Type_t,3>( op( lhs[0], rhs[0] ) ,
                                    op( lhs[1], rhs[1] ) ,
                                    op( lhs[2], rhs[2] ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< AntiSymTensor<T1,3> , T2 , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static AntiSymTensor<Type_t,3>
  apply(const AntiSymTensor<T1,3>& lhs, T2 rhs, OP op)
  {
    return AntiSymTensor<Type_t,3>( op( lhs[0], rhs ) ,
                                    op( lhs[1], rhs ) ,
                                    op( lhs[2], rhs ) );
  }
};

template<class T1, class T2, class OP>
struct OTBinary< T1, AntiSymTensor<T2,3> , OP >
{
  typedef typename BinaryReturn<T1,T2,OP>::Type_t Type_t;
  inline static AntiSymTensor<Type_t,3>
  apply(T1 lhs, const AntiSymTensor<T2,3>& rhs, OP op)
  {
    return AntiSymTensor<Type_t,3>( op( lhs, rhs[0] ) ,
                                    op( lhs, rhs[1] ) ,
                                    op( lhs, rhs[2] ) );
  }
};


//////////////////////////////////////////////////////
//
// determinant: generalized
//
//////////////////////////////////////////////////////
template <class T, unsigned D>
inline typename Tensor<T,D>::Type_t det(const Tensor<T,D>& a)
{
  // to implement the general case here
  return  0;
}

//////////////////////////////////////////////////////
// specialized for D=1
//////////////////////////////////////////////////////
template <class T>
inline typename Tensor<T,1>::Type_t det(const Tensor<T,1>& a)
{
  return  a(0,0);
}

//////////////////////////////////////////////////////
// specialized for D=2
//////////////////////////////////////////////////////
template <class T>
inline typename Tensor<T,2>::Type_t det(const Tensor<T,2>& a)
{
  return  a(0,0)*a(1,1)-a(0,1)*a(1,0);
}

//////////////////////////////////////////////////////
// specialized for D=3
//////////////////////////////////////////////////////
template <class T>
inline typename Tensor<T,3>::Type_t det(const Tensor<T,3>& a)
{
  return  a(0,0)*(a(1,1)*a(2,2)-a(1,2)*a(2,1))
          +a(0,1)*(a(1,2)*a(2,0)-a(1,0)*a(2,2))
          +a(0,2)*(a(1,0)*a(2,1)-a(1,1)*a(2,0));
}

//////////////////////////////////////////////////////
//
// inverse: generalized
// A*B = I, * being the matrix multiplication  I(i,j) = sum_k A(i,k)*B(k,j)
//
//////////////////////////////////////////////////////
template <class T, unsigned D>
inline Tensor<T,D> inverse(const Tensor<T,D>& a)
{
  return Tensor<T,D>();
}

//////////////////////////////////////////////////////
// specialized for D=1
//////////////////////////////////////////////////////
template <class T>
inline Tensor<T,1> inverse(const Tensor<T,1>& a)
{
  return Tensor<T,1>(1.0/a(0,0));
}

//////////////////////////////////////////////////////
// specialized for D=2
//////////////////////////////////////////////////////
template <class T>
inline Tensor<T,2> inverse(const Tensor<T,2>& a)
{
  T vinv=1/det(a);
  return Tensor<T,2>(vinv*a(1,1), -vinv*a(0,1), -vinv*a(1,0), vinv*a(0,0));
}

//////////////////////////////////////////////////////
// specialized for D=3
//////////////////////////////////////////////////////
template <class T>
inline Tensor<T,3> inverse(const Tensor<T,3>& a)
{
  T vinv=1/det(a);
//   return Tensor<T,3>(vinv*(a(1,1)*a(2,2)-a(1,2)*a(2,1)),
//                      vinv*(a(1,2)*a(2,0)-a(1,0)*a(2,2)),
//                      vinv*(a(1,0)*a(2,1)-a(1,1)*a(2,0)),
//                      vinv*(a(2,1)*a(0,2)-a(2,2)*a(0,1)),
//                      vinv*(a(2,2)*a(0,0)-a(2,0)*a(0,2)),
//                      vinv*(a(2,0)*a(0,1)-a(2,1)*a(0,0)),
//                      vinv*(a(0,1)*a(1,2)-a(0,2)*a(1,1)),
//                      vinv*(a(0,2)*a(1,0)-a(0,0)*a(1,2)),
//                      vinv*(a(0,0)*a(1,1)-a(0,1)*a(1,0)));
  return Tensor<T,3>(vinv*(a(1,1)*a(2,2)-a(1,2)*a(2,1)),
                     vinv*(a(2,1)*a(0,2)-a(2,2)*a(0,1)),
                     vinv*(a(0,1)*a(1,2)-a(0,2)*a(1,1)),
                     vinv*(a(1,2)*a(2,0)-a(1,0)*a(2,2)),
                     vinv*(a(2,2)*a(0,0)-a(2,0)*a(0,2)),
                     vinv*(a(0,2)*a(1,0)-a(0,0)*a(1,2)),
                     vinv*(a(1,0)*a(2,1)-a(1,1)*a(2,0)),
                     vinv*(a(2,0)*a(0,1)-a(2,1)*a(0,0)),
                     vinv*(a(0,0)*a(1,1)-a(0,1)*a(1,0)));
//   int i,j,i1,i2,j1,j2;
//   int cyclic[]={1,2,0};
//   for(i=0;i<3;i++){
//     i1 = cyclic[i];
//     i2 = cyclic[i1];
//     for(j=0;j<3;j++){
//       j1 = cyclic[j];
//       j2 = cyclic[j1];
//       x(i,j) = a(i1,j1)*a(i2,j2) - a(i1,j2)*a(i2,j1);
//     }
//  }
//  typename T_t::Element_t
//  detinv= 1.0e0/(a(0,0)*x(0,0)+a(1,0)*x(1,0)+a(2,0)*x(2,0));
//  return detinv*x;
}

//////////////////////////////////////////////////////////////////////
//
// Specializations for Tensor dot Tensor
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, unsigned D>
struct OTDot< Tensor<T1,D> , Tensor<T2,D> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Tensor<Type_t,D>
  apply(const Tensor<T1,D>& lhs, const Tensor<T2,D>& rhs)
  {
    Tensor<Type_t,D> res = Tensor<Type_t,D>::DontInitialize();
    for (unsigned int i=0; i<D; ++i)
      for (unsigned int j=0; j<D; ++j)
      {
        Type_t sum = lhs(i,0) * rhs(0,j);
        for (unsigned int k=1; k<D; ++k)
          sum += lhs(i,k) * rhs(k,j);
        res(i,j) = sum;
      }
    return res;
  }
};

template<class T1, class T2>
struct OTDot< Tensor<T1,1> , Tensor<T2,1> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Tensor<Type_t,1>
  apply(const Tensor<T1,1>& lhs, const Tensor<T2,1>& rhs)
  {
    return Tensor<Type_t,1>(lhs[0]*rhs[0]);
  }
};

template<class T1, class T2>
struct OTDot< Tensor<T1,2> , Tensor<T2,2> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Tensor<Type_t,2>
  apply(const Tensor<T1,2>& lhs, const Tensor<T2,2>& rhs)
  {
    return Tensor<Type_t,2>(lhs(0,0)*rhs(0,0) + lhs(0,1)*rhs(1,0),
                            lhs(0,0)*rhs(0,1) + lhs(0,1)*rhs(1,1),
                            lhs(1,0)*rhs(0,0) + lhs(1,1)*rhs(1,0),
                            lhs(1,0)*rhs(0,1) + lhs(1,1)*rhs(1,1));
  }
};

template<class T1, class T2>
struct OTDot< Tensor<T1,3> , Tensor<T2,3> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Tensor<Type_t,3>
  apply(const Tensor<T1,3>& lhs, const Tensor<T2,3>& rhs)
  {
    return Tensor<Type_t,3>( lhs(0,0)*rhs(0,0) + lhs(0,1)*rhs(1,0) + lhs(0,2)*rhs(2,0) ,
                             lhs(0,0)*rhs(0,1) + lhs(0,1)*rhs(1,1) + lhs(0,2)*rhs(2,1) ,
                             lhs(0,0)*rhs(0,2) + lhs(0,1)*rhs(1,2) + lhs(0,2)*rhs(2,2) ,
                             lhs(1,0)*rhs(0,0) + lhs(1,1)*rhs(1,0) + lhs(1,2)*rhs(2,0) ,
                             lhs(1,0)*rhs(0,1) + lhs(1,1)*rhs(1,1) + lhs(1,2)*rhs(2,1) ,
                             lhs(1,0)*rhs(0,2) + lhs(1,1)*rhs(1,2) + lhs(1,2)*rhs(2,2) ,
                             lhs(2,0)*rhs(0,0) + lhs(2,1)*rhs(1,0) + lhs(2,2)*rhs(2,0) ,
                             lhs(2,0)*rhs(0,1) + lhs(2,1)*rhs(1,1) + lhs(2,2)*rhs(2,1) ,
                             lhs(2,0)*rhs(0,2) + lhs(2,1)*rhs(1,2) + lhs(2,2)*rhs(2,2) );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations for SymTensor dot SymTensor
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, unsigned D>
struct OTDot< SymTensor<T1,D> , SymTensor<T2,D> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Tensor<Type_t,D>
  apply(const SymTensor<T1,D>& lhs, const SymTensor<T2,D>& rhs)
  {
    Tensor<Type_t,D> res = Tensor<Type_t,D>::DontInitialize();
    for (unsigned int i=0; i<D; ++i)
      for (unsigned int j=i; j<D; ++j)
      {
        Type_t sum = lhs.HL(i,0) * rhs.HL(j,0);
        int k=1;
        for ( ; k<i; ++k )
          sum += lhs.HL(i,k) * rhs.HL(j,k);
        for ( ; k<j; ++k )
          sum += lhs.HL(k,i) * rhs.HL(j,k);
        for ( ; k<D; ++k )
          sum += lhs.HL(k,i) * rhs.HL(k,j);
        res(i,j) = sum;
      }
    return res;
  }
};

template<class T1, class T2>
struct OTDot< SymTensor<T1,1> , SymTensor<T2,1> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Tensor<Type_t,1>
  apply(const SymTensor<T1,1>& lhs, const SymTensor<T2,1>& rhs)
  {
    return Tensor<Type_t,1>(lhs[0]*rhs[0]);
  }
};

template<class T1, class T2>
struct OTDot< SymTensor<T1,2> , SymTensor<T2,2> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Tensor<Type_t,2>
  apply(const SymTensor<T1,2>& lhs, const SymTensor<T2,2>& rhs)
  {
    return Tensor<Type_t,2>(lhs(0,0)*rhs(0,0) + lhs(0,1)*rhs(1,0),
                            lhs(0,0)*rhs(0,1) + lhs(0,1)*rhs(1,1),
                            lhs(1,0)*rhs(0,0) + lhs(1,1)*rhs(1,0),
                            lhs(1,0)*rhs(0,1) + lhs(1,1)*rhs(1,1));
  }
};

template<class T1, class T2>
struct OTDot< SymTensor<T1,3> , SymTensor<T2,3> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Tensor<Type_t,3>
  apply(const SymTensor<T1,3>& lhs, const SymTensor<T2,3>& rhs)
  {
    return Tensor<Type_t,3>( lhs(0,0)*rhs(0,0) + lhs(0,1)*rhs(1,0) + lhs(0,2)*rhs(2,0) ,
                             lhs(0,0)*rhs(0,1) + lhs(0,1)*rhs(1,1) + lhs(0,2)*rhs(2,1) ,
                             lhs(0,0)*rhs(0,2) + lhs(0,1)*rhs(1,2) + lhs(0,2)*rhs(2,2) ,
                             lhs(1,0)*rhs(0,0) + lhs(1,1)*rhs(1,0) + lhs(1,2)*rhs(2,0) ,
                             lhs(1,0)*rhs(0,1) + lhs(1,1)*rhs(1,1) + lhs(1,2)*rhs(2,1) ,
                             lhs(1,0)*rhs(0,2) + lhs(1,1)*rhs(1,2) + lhs(1,2)*rhs(2,2) ,
                             lhs(2,0)*rhs(0,0) + lhs(2,1)*rhs(1,0) + lhs(2,2)*rhs(2,0) ,
                             lhs(2,0)*rhs(0,1) + lhs(2,1)*rhs(1,1) + lhs(2,2)*rhs(2,1) ,
                             lhs(2,0)*rhs(0,2) + lhs(2,1)*rhs(1,2) + lhs(2,2)*rhs(2,2) );
  }
};
}

#endif // OHMMS_TENSOR_OPERATORS_H


