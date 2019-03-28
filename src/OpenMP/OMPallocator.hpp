//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file OMPallocator.hpp
 */
#ifndef QMCPLUSPLUS_OPENMP_ALLOCATOR_H
#define QMCPLUSPLUS_OPENMP_ALLOCATOR_H

#include <memory>
#include "config.h"

namespace qmcplusplus
{
  template<typename T, class HostAllocator=std::allocator<T>>
  struct OMPallocator: public HostAllocator
  {
    using value_type    = typename HostAllocator::value_type;
    using size_type     = typename HostAllocator::size_type;
    using pointer       = typename HostAllocator::pointer;
    using const_pointer = typename HostAllocator::const_pointer;

    OMPallocator() = default;
    template <class U, class V> OMPallocator(const OMPallocator<U, V>&) {}
    template <class U, class V> struct rebind { typedef OMPallocator<U, V> other; };

    value_type* allocate(std::size_t n)
    {
      value_type* pt = HostAllocator::allocate(n);
      PRAGMA_OFFLOAD("omp target enter data map(alloc:pt[0:n])")
      return pt;
    }

    void deallocate(value_type* pt, std::size_t n)
    {
      PRAGMA_OFFLOAD("omp target exit data map(delete:pt[0:n])")
      HostAllocator::deallocate(pt,n);
    }
  };
}

#endif
