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
    template <class U> OMPallocator(const OMPallocator<U>&) {}
    template <class U> struct rebind { typedef OMPallocator<U> other; };

    value_type* allocate(std::size_t n, int device_id = 0)
    {
      #pragma omp target enter data map(alloc:pt[0:n]) device(device_id)
      value_type* pt = HostAllocator::allocate(n);
      return pt;
    }

    void deallocate(value_type* pt, std::size_t, int device_id = 0) {
      #pragma omp target exit data map(delete:pt) device(device_id)
      free(pt);
    }
  };

  template <class HOST1, class HOST2>
  bool operator==(const OMPallocator<HOST1>&, const OMPallocator<HOST2>&) { return HOST1==HOST2; }
  template <class HOST1, class HOST2>
  bool operator!=(const OMPallocator<HOST1>&, const OMPallocator<HOST2>&) { return HOST1!=HOST2; }
}

#endif
