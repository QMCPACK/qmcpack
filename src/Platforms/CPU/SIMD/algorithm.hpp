//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by:
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file algorithm.hpp
 *
 * SIMD version of functions in algorithm
 */
#ifndef QMCPLUSPLUS_SIMD_ALGORITHM_HPP
#define QMCPLUSPLUS_SIMD_ALGORITHM_HPP

#include <complex>

namespace qmcplusplus
{
namespace simd
{
template<typename T1, typename T2>
inline T2 accumulate_n(const T1* restrict in, size_t n, T2 res)
{
#pragma omp simd reduction(+ : res)
  for (int i = 0; i < n; ++i)
    res += in[i];
  return res;
}

/** copy function using memcpy
     * @param target starting address of the target
     * @param source starting address of the source
     * @param n size of the data to copy
     */
template<typename T1, typename T2>
inline void copy(T1* restrict target, const T2* restrict source, size_t n)
{
  for (size_t i = 0; i < n; ++i)
    target[i] = static_cast<T1>(source[i]);
}

/** copy function using memcpy
     * @param target starting address of the target
     * @param source starting address of the source
     * @param n size of the data to copy
     */
template<typename T>
inline void copy(T* restrict target, const T* restrict source, size_t n)
{
  memcpy(target, source, sizeof(T) * n);
}

/** copy complex to two real containers */
template<typename T1, typename T2>
inline void copy(T1* restrict target_r, T1* restrict target_i, const std::complex<T2>* restrict source, size_t n)
{
  for (int i = 0; i < n; ++i)
  {
    *target_r++ = static_cast<T1>(source[i].real());
    *target_i++ = static_cast<T1>(source[i].imag());
  }
}

template<typename T>
inline void accumulate_phases(const int& n,
                              const std::complex<T>* restrict x,
                              const std::complex<T>* restrict y,
                              T& rN,
                              T& iN,
                              T& riN)
{
  for (int i = 0; i < n; ++i)
  {
    T tr = x[i].real() * y[i].real() - x[i].imag() * y[i].imag();
    T ti = x[i].real() * y[i].imag() + x[i].imag() * y[i].real();
    rN += tr * tr;
    iN += ti * ti;
    riN += tr * ti;
  } //
}

/** transpose of A(m,n) to B(n,m)
     * @param A starting address, A(m,lda)
     * @param m number of A rows
     * @param lda stride of A's row
     * @param B starting address B(n,ldb)
     * @param n number of B rows
     * @param ldb stride of B's row
     *
     * Blas-like interface
     */
template<typename T, typename TO>
inline void transpose(const T* restrict A, size_t m, size_t lda, TO* restrict B, size_t n, size_t ldb)
{
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < m; ++j)
      B[i * ldb + j] = A[j * lda + i];
}

/** copy of A(m,n) to B(m,n)
     * @param A starting address, A(m,lda)
     * @param m number of A rows
     * @param lda stride of A's row
     * @param B starting address B(m,ldb)
     * @param m number of B rows
     * @param ldb stride of B's row
     *
     * Blas-like interface
     */
template<typename T, typename TO>
inline void remapCopy(size_t m, size_t n, const T* restrict A, size_t lda, TO* restrict B, size_t ldb)
{
  for (size_t j = 0; j < m; ++j)
    for (size_t i = 0; i < n; ++i)
      B[j * ldb + i] = A[j * lda + i];
}
} // namespace simd
} // namespace qmcplusplus
#endif
