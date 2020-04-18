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
#include <type_traits>
#include "config.h"

namespace qmcplusplus
{
template<typename T, class HostAllocator = std::allocator<T>>
struct OMPallocator : public HostAllocator
{
  using value_type    = typename HostAllocator::value_type;
  using size_type     = typename HostAllocator::size_type;
  using pointer       = typename HostAllocator::pointer;
  using const_pointer = typename HostAllocator::const_pointer;

  OMPallocator() = default;
  template<class U, class V>
  OMPallocator(const OMPallocator<U, V>&)
  {}
  template<class U, class V>
  struct rebind
  {
    typedef OMPallocator<U, V> other;
  };

  value_type* allocate(std::size_t n)
  {
    static_assert(std::is_same<T, value_type>::value, "OMPallocator and HostAllocator data types must agree!");
    value_type* pt = HostAllocator::allocate(n);
    PRAGMA_OFFLOAD("omp target enter data map(alloc:pt[0:n])")
    return pt;
  }

  void deallocate(value_type* pt, std::size_t n)
  {
    PRAGMA_OFFLOAD("omp target exit data map(delete:pt[0:n])")
    HostAllocator::deallocate(pt, n);
  }
};

template<typename T>
auto* getContainerDevicePtr(T& dataset)
{
  auto* host_ptr = dataset.data();
  typename T::value_type* device_ptr;
  PRAGMA_OFFLOAD("omp target data use_device_ptr(host_ptr)")
  {
    device_ptr = host_ptr;
  }
  return device_ptr;
}

} // namespace qmcplusplus

#endif
