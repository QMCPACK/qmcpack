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

#include <config.h>
#include <vector>
#include <cstdlib>

#if (__cplusplus >= 201103L)
#if defined(__INTEL_COMPILER)
 #include <tbb/cache_aligned_allocator.h>
#else
 #include "simd/Mallocator.hpp"
#endif

namespace qmcplusplus
{
  template<class T>
#if defined(__INTEL_COMPILER)
    using aligned_allocator=tbb::cache_aligned_allocator<T>;
#elif __bgq__
    using aligned_allocator=std::allocator<T>;
#else
    using aligned_allocator=qmcplusplus::Mallocator<T, QMC_CLINE>;
#endif
  template<class T>
    using aligned_vector = std::vector<T,aligned_allocator<T> >;

}
#else
namespace qmcplusplus
{
  /** dummy inherited class to use aligned_vector<T> */
  template<class T> struct aligned_vector: public std::vector<T> { };

  /** dummy inherited class to use aligned_allocator<T> */
  template<class T> struct aligned_allocator: public std::allocator<T> { };
}
#endif

template<typename T> inline size_t getAlignedSize(size_t n)
{
  CONSTEXPR size_t ND=QMC_CLINE/sizeof(T);
  return ((n+ND-1)/ND)*ND;
}


#endif
