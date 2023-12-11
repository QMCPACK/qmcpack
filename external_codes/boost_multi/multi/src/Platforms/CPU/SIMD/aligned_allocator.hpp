//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file allocator.hpp
 */
#ifndef QMCPLUSPLUS_ALLOCATOR_H
#define QMCPLUSPLUS_ALLOCATOR_H

#include <vector>
#include <cstdlib>
#include "config.h"
#include "Mallocator.hpp"

#if defined(__INTEL_COMPILER)
  #define ASSUME_ALIGNED(x) __assume_aligned(x,QMC_SIMD_ALIGNMENT)
#elif defined(__GNUC__) && !defined(__ibmxl__)
  #define ASSUME_ALIGNED(x) (x) = (__typeof__(x)) __builtin_assume_aligned(x,QMC_SIMD_ALIGNMENT)
#else
  #define ASSUME_ALIGNED(x)
#endif

namespace qmcplusplus
{
template<class T, size_t ALIGN = QMC_SIMD_ALIGNMENT>
using aligned_allocator = Mallocator<T, ALIGN>;
template<class T>
using aligned_vector = std::vector<T, aligned_allocator<T>>;

} // namespace qmcplusplus

/** return size in T's of allocated aligned memory
 */
template<typename T, size_t ALIGN = QMC_SIMD_ALIGNMENT>
inline size_t getAlignedSize(size_t n)
{
  constexpr size_t ND = ALIGN / sizeof(T);
  static_assert(ALIGN % sizeof(T) == 0, "getAlignedSize ALIGN must be a multiple of sizeof(T)");
  return ((n + ND - 1) / ND) * ND;
}

template<typename T, size_t ALIGN = QMC_SIMD_ALIGNMENT>
inline size_t getAlignment()
{
  static_assert(ALIGN % sizeof(T) == 0, "getAlignedSize ALIGN must be a multiple of sizeof(T)");
  return ALIGN / sizeof(T);
}

#endif
