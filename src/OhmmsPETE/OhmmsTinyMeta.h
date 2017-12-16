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


#ifndef OHMMS_TINY_META_H
#define OHMMS_TINY_META_H


namespace qmcplusplus
{
/* \note
 *  For optimization, operators for TinyVector and Tensor classes are specalized for D.
 *  Instead of using PETE generated operators, the operators for TinyVector, Tensor and TinyMatrix
 *  classes are explicitly implemented in several files.
 *  The collections are rather arbitrary but follow the r1/src/AppTypes pattern closely.
 *  This file is included by TinyVector.h, Tensor.h, and TinyMatrix.h.
 *  Therefore, users do not have to include this file.
 */

// foward declarations
template<class T, unsigned D> class TinyVector;
template<class T, unsigned D> class Tensor;
template<class T, unsigned D> class SymTensor;
template<class T, unsigned D> class AntiSymTensor;
template<class T, unsigned D1, unsigned D2> class TinyMatrix;

// generic OTAssign and OTBinary functors (OT = OhmmsTiny)
template<class T1, class T2, class OP> struct OTAssign {};
template<class T1, class T2, class OP> struct OTBinary {};

// generic Dot-product functors
template<class T1, class T2> struct OTDot {};
}

////////////////////////////////////////////////////////////////////////////////
// TinyVectorOps.h        assignment/unary and binary operators  for TinyVector
// TensorOps.h            assignment/unary and binary operators  for Tensors
// TinyVectorTensorOps.h  binary operators for TinyVector and Tensor combinations
// TinyMatrixOps.h        assignment/unary and binary operators  for TinyMatrix
////////////////////////////////////////////////////////////////////////////////
#define PAssert
#include "OhmmsPETE/TinyVectorOps.h"
#include "OhmmsPETE/TensorOps.h"
#include "OhmmsPETE/TinyVectorTensorOps.h"
#include "OhmmsPETE/TinyMatrixOps.h"

// macros to generate a set of binary and unary combintions for each operator
// e.g., OHMMS_META_BINARY_OPERATORS(TinyVector,operator+,OpAdd)
// generates
// a) TinyVector + TinyVector
// b) TinyVector + scalar
// c) scalar     + TinyVector
// TinyVector.h and Tensor.h choose appropriate operators.

#define OHMMS_META_BINARY_OPERATORS(TENT,FUNC,TAG)                            \
                                                                              \
template <class T1, class T2, unsigned D>                                     \
struct BinaryReturn< TENT<T1,D>, TENT<T2,D>, TAG >                            \
{                                                                             \
  typedef TENT<typename BinaryReturn<T1,T2,TAG>::Type_t, D> Type_t;           \
};                                                                            \
                                                                              \
template <class T1, class T2, unsigned D>                                     \
inline                                                                        \
typename BinaryReturn< TENT<T1,D>, TENT<T2,D>, TAG >::Type_t                  \
FUNC( const TENT<T1,D>& v1, const TENT<T2,D>& v2 )                            \
{                                                                             \
  return OTBinary<TENT<T1,D>,TENT<T2,D>,TAG>::apply(v1,v2,TAG());             \
}                                                                             \
                                                                              \
template <class T1, class T2, unsigned D>                                     \
struct BinaryReturn< TENT<T1,D>, T2, TAG >                                    \
{                                                                             \
  typedef TENT< typename BinaryReturn<T1,T2,TAG>::Type_t, D > Type_t;         \
};                                                                            \
                                                                              \
template <class T1, class T2, unsigned D>                                     \
struct BinaryReturn< T1, TENT<T2,D>, TAG >                                    \
{                                                                             \
  typedef TENT< typename BinaryReturn<T1,T2,TAG>::Type_t, D > Type_t;         \
};                                                                            \
                                                                              \
template <class T1, class T2, unsigned D>                                     \
inline                                                                        \
typename BinaryReturn< TENT<T1,D>, T2, TAG >::Type_t                          \
FUNC( const TENT<T1,D>& v1, const T2& x )                                     \
{                                                                             \
  return OTBinary<TENT<T1,D>,T2,TAG>::apply(v1,x,TAG());                      \
}                                                                             \
                                                                              \
template <class T1, class T2, unsigned D>                                     \
inline                                                                        \
typename BinaryReturn< T1, TENT<T2,D>, TAG >::Type_t                          \
FUNC( const T1& x, const TENT<T2,D>& v2)                                      \
{                                                                             \
  return OTBinary<T1,TENT<T2,D>,TAG>::apply(x,v2,TAG());                      \
}                                                                             \
 
#define OHMMS_META_ACCUM_OPERATORS(TENT,FUNC,TAG)                             \
                                                                              \
template <class T1, class T2, unsigned D>                                     \
inline TENT<T1,D>&                                                            \
FUNC( TENT<T1,D>& v1, const TENT<T2,D>& v2 )                                  \
{                                                                             \
  OTAssign<TENT<T1,D>,TENT<T2,D>,TAG>::apply(v1,v2,TAG());                    \
  return v1;                                                                  \
}                                                                             \
                                                                              \
template <class T1, class T2, unsigned D>                                     \
inline TENT<T1,D>&                                                            \
FUNC( TENT<T1,D>& v1, const T2& v2 )                                          \
{                                                                             \
  OTAssign<TENT<T1,D>,T2,TAG>::apply(v1,v2,TAG());                            \
  return v1;                                                                  \
}                                                                             \
 
#define OHMMS_TINYMAT_BINARY_OPERATORS(TENT,FUNC,TAG)                         \
                                                                              \
template <class T1, class T2, unsigned D1, unsigned D2>                       \
struct BinaryReturn< TENT<T1,D1,D2>, TENT<T2,D1,D2>, TAG >                    \
{                                                                             \
  typedef TENT<typename BinaryReturn<T1,T2,TAG>::Type_t,D1,D2> Type_t;        \
};                                                                            \
                                                                              \
template <class T1, class T2, unsigned D1, unsigned D2>                       \
inline                                                                        \
typename BinaryReturn< TENT<T1,D1,D2>, TENT<T2,D1,D2>, TAG >::Type_t          \
FUNC( const TENT<T1,D1,D2>& v1, const TENT<T2,D1,D2>& v2 )                    \
{                                                                             \
  return OTBinary<TENT<T1,D1,D2>,TENT<T2,D1,D2>,TAG>::apply(v1,v2,TAG());     \
}                                                                             \
                                                                              \
template <class T1, class T2, unsigned D1, unsigned D2>                       \
struct BinaryReturn< TENT<T1,D1,D2>, T2, TAG >                                \
{                                                                             \
  typedef TENT< typename BinaryReturn<T1,T2,TAG>::Type_t, D1,D2> Type_t;      \
};                                                                            \
                                                                              \
template <class T1, class T2, unsigned D1, unsigned D2>                       \
struct BinaryReturn< T1, TENT<T2,D1,D2>, TAG >                                \
{                                                                             \
  typedef TENT< typename BinaryReturn<T1,T2,TAG>::Type_t, D1,D2> Type_t;      \
};                                                                            \
                                                                              \
template <class T1, class T2, unsigned D1, unsigned D2>                       \
inline                                                                        \
typename BinaryReturn< TENT<T1,D1,D2>, T2, TAG >::Type_t                      \
FUNC( const TENT<T1,D1,D2>& v1, const T2& x )                                 \
{                                                                             \
  return OTBinary<TENT<T1,D1,D2>,T2,TAG>::apply(v1,x,TAG());                  \
}                                                                             \
                                                                              \
template <class T1, class T2, unsigned D1, unsigned D2>                       \
inline                                                                        \
typename BinaryReturn< T1, TENT<T2,D1,D2>, TAG >::Type_t                      \
FUNC( const T1& x, const TENT<T2,D1,D2>& v2)                                  \
{                                                                             \
  return OTBinary<T1,TENT<T2,D1,D2>,TAG>::apply(x,v2,TAG());                  \
}                                                                             \
 
#define OHMMS_TINYMAT_ACCUM_OPERATORS(TENT,FUNC,TAG)                          \
                                                                              \
template <class T1, class T2, unsigned D1, unsigned D2>                       \
inline TENT<T1, D1, D2>&                                                      \
FUNC( TENT<T1, D1, D2>& v1, const TENT<T2, D1, D2>& v2 )                      \
{                                                                             \
  OTAssign<TENT<T1,D1,D2>,TENT<T2,D1,D2>,TAG>::apply(v1,v2,TAG());            \
  return v1;                                                                  \
}                                                                             \
                                                                              \
template <class T1, class T2, unsigned D1, unsigned D2>                       \
inline TENT<T1,D1,D2>&                                                        \
FUNC( TENT<T1,D1,D2>& v1, const T2& v2 )                                      \
{                                                                             \
  OTAssign<TENT<T1,D1,D2>,T2,TAG>::apply(v1,v2,TAG());                        \
  return v1;                                                                  \
}                                                                             \
 
#endif // OHMMS_TINY_META_H

