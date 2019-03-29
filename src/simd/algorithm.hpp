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

///inner product
template<typename T1, typename T2, typename T3>
inline T3 inner_product_n(const T1* restrict a, const T2* restrict b, int n, T3 res)
{
  for (int i = 0; i < n; ++i)
    res += a[i] * b[i];
  return res;
}

} // namespace simd
} // namespace qmcplusplus
#endif
