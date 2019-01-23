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
#include <stdexcept>

namespace qmcplusplus
{
  template<typename T, size_t ALIGN>
  struct Mallocator
  {
    typedef T         value_type;
    typedef size_t    size_type;
    typedef T*        pointer;
    typedef const T*  const_pointer;

    Mallocator() = default;
    template <class U> Mallocator(const Mallocator<U,ALIGN>&) {}

    template <class U> struct rebind { typedef Mallocator<U, ALIGN> other; };

    T* allocate(std::size_t n)
    {
      void* pt(nullptr);
#if __GLIBC__ == 2 && __GLIBC_MINOR__ >= 16
      std::size_t asize = n * sizeof(T);
      std::size_t amod = asize % ALIGN;
      if (amod != 0) asize += ALIGN - amod;
      // as per C++11 standard asize must be an integral multiple of ALIGN
      // or behavior is undefined.  Some implementation support all positive
      // values of asize but the standard has not been amended
      // This is also not guaranteed threadsafe until C++17
      pt = aligned_alloc(ALIGN,asize);
#else
      posix_memalign(&pt, ALIGN, n*sizeof(T));
#endif
      if ( pt == nullptr )
        throw std::runtime_error("Allocation failed in Mallocator, requested size in bytes = " + std::to_string(n*sizeof(T)));
      return static_cast<T*>(pt);
    }

    void deallocate(T* p, std::size_t) {
      free(p);
    }
  };

  template <class T1, size_t ALIGN1, class T2, size_t ALIGN2>
  bool operator==(const Mallocator<T1,ALIGN1>&, const Mallocator<T2,ALIGN2>&) { return ALIGN1==ALIGN2; }
  template <class T1, size_t ALIGN1, class T2, size_t ALIGN2>
  bool operator!=(const Mallocator<T1,ALIGN1>&, const Mallocator<T2,ALIGN2>&) { return ALIGN1!=ALIGN2; }
}

#endif
