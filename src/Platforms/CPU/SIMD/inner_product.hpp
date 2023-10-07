//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file inner_product.h
 *
 * Inline dot and gemv functions. 
 * On modern processors with good compilers, there is no need to use blas1, e.g., ddot, but 
 * strange things can happen on some machines.
 * SIMD specialization can be implemented here.
 */
#ifndef QMCPLUSPLUS_INNER_PRODUCT_HPP
#define QMCPLUSPLUS_INNER_PRODUCT_HPP

#include "OhmmsPETE/TinyVector.h"

namespace qmcplusplus
{
namespace simd
{
///inner product
template<typename T1, typename T2, typename T3>
inline T3 inner_product_n(const T1* restrict a, const T2* restrict b, int n, T3 res)
{
  for (int i = 0; i < n; ++i)
    res += a[i] * b[i];
  return res;
}

/** dot product
     * @param a starting address of an array of type T
     * @param b starting address of an array of type T
     * @param n size
     * @param res initial value, default=0.0
     * @return  \f$ res = \sum_i a[i] b[i]\f$
     *
     * same as inner_product(a,a+n,b,0.0)
     * The arguments of this inline function are different from BLAS::dot
     * This is more efficient than BLAS::dot due to the function overhead,
     * if a compiler knows how to inline.
     */
template<typename T>
inline T dot(const T* restrict a, const T* restrict b, int n, T res = T())
{
  for (int i = 0; i < n; i++)
    res += a[i] * b[i];
  return res;
}

/** inline dot product 
     * @param a starting address of an array of type T
     * @param b starting address of an array of type TinyVector<T,D>
     * @param n size
     * @return \f$ {\bf v} = \sum_i a[i] {\bf b}[i]\f$
     */
template<class T, unsigned D>
inline Tensor<T, D> dot(const T* a, const Tensor<T, D>* b, int n)
{
  Tensor<T, D> res;
  //#if defined(USE_GEMV_FOR_G_DOT_V)
  //        BLAS::gemv(n,3,b->data(),a,res.data());
  //#else
  for (int i = 0; i < n; i++)
    res += a[i] * b[i];
  //#endif
  return res;
}

template<class T, unsigned D>
inline TinyVector<T, D> dot(const T* a, const TinyVector<T, D>* b, int n)
{
  TinyVector<T, D> res;
#if defined(USE_GEMV_FOR_G_DOT_V)
  BLAS::gemv(n, 3, b->data(), a, res.data());
#else
  for (int i = 0; i < n; i++)
    res += a[i] * b[i];
#endif
  return res;
}

/** x*y dot product of two vectors using the same argument list for blas::dot
     * @param n size
     * @param x starting address of x 
     * @param incx stride of x
     * @param y starting address of y 
     * @param incx stride of y
     * @param return \f$\sum_i x[i+=incx]*y[i+=incy]\f$
     */
template<typename T>
inline T dot(int n, const T* restrict x, int incx, const T* restrict y, int incy)
{
  const int xmax = incx * n;
  T res          = 0.0;
  for (int ic = 0, jc = 0; ic < xmax; ic += incx, jc += incy)
    res += x[ic] * y[jc];
  return res;
}

} // namespace simd
} // namespace qmcplusplus
#endif
