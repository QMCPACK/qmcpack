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

namespace qmcplusplus
{
  template<typename T, size_t Align>
  struct Mallocator
  {
    typedef T value_type;
    typedef size_t size_type;
    Mallocator() = default;
    template <class U> Mallocator(const Mallocator<U,Align>&) {}

    template <class U> struct rebind { typedef Mallocator<U, Align> other; };

    T* allocate(std::size_t n) {
#if __GLIBC__ == 2 && __GLIBC_MINOR__ >= 16
      return static_cast<T*>(aligned_alloc(Align,n*sizeof(T)));
#else
      void* pt;
      posix_memalign(&pt, Align, n*sizeof(T));
      return static_cast<T*>(pt);
#endif
    }
    void deallocate(T* p, std::size_t) { free(p); }
  };
}

//  template<class T>
//   using aligned_allocator=Mallocator<T,QMC_CLINE>;
#endif
