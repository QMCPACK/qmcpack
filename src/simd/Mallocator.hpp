//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file Mallocator.hpp
 */
#ifndef QMCPLUSPLUS_ALIGNED_ALLOCATOR_H
#define QMCPLUSPLUS_ALIGNED_ALLOCATOR_H

#include <cstdlib>
#include <string>
#include <stdexcept>
#include <string>

namespace qmcplusplus
{
template<typename T, size_t ALIGN>
struct Mallocator
{
  typedef T value_type;
  typedef size_t size_type;
  typedef T* pointer;
  typedef const T* const_pointer;

  Mallocator() = default;
  template<class U>
  Mallocator(const Mallocator<U, ALIGN>&)
  {}

  template<class U>
  struct rebind
  {
    typedef Mallocator<U, ALIGN> other;
  };

  T* allocate(std::size_t n)
  {
    void* pt(nullptr);
    std::size_t asize = n * sizeof(T);
    std::size_t amod  = asize % ALIGN;
    if (amod != 0)
      asize += ALIGN - amod;

#if __STDC_VERSION__ >= 201112L || (__GLIBC__ == 2 && __GLIBC_MINOR__ >= 16)
    // as per C11 standard asize must be an integral multiple of ALIGN
    // or behavior is undefined.  Some implementations support all positive
    // values of asize but the standard has not been amended
    // This is also not guaranteed threadsafe until it appeared in
    // the C++17 standard.
    pt = aligned_alloc(ALIGN, asize);
#else
    // While posix memalign can deal with asize violating the C11 standard
    // assumptions made later by our simd code namely copyn require allocation
    // of the entire aligned block to avoid heap buffer read overflows later
    posix_memalign(&pt, ALIGN, asize);
#endif
    if (pt == nullptr)
      throw std::runtime_error("Allocation failed in Mallocator, requested size in bytes = " +
                               std::to_string(n * sizeof(T)));
    return static_cast<T*>(pt);
  }

  void deallocate(T* p, std::size_t) { free(p); }
};

template<class T1, size_t ALIGN1, class T2, size_t ALIGN2>
bool operator==(const Mallocator<T1, ALIGN1>&, const Mallocator<T2, ALIGN2>&)
{
  return ALIGN1 == ALIGN2;
}
template<class T1, size_t ALIGN1, class T2, size_t ALIGN2>
bool operator!=(const Mallocator<T1, ALIGN1>&, const Mallocator<T2, ALIGN2>&)
{
  return ALIGN1 != ALIGN2;
}
} // namespace qmcplusplus

#endif
