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
#include <atomic>
#include "config.h"
#include "allocator_traits.hpp"

namespace qmcplusplus
{

extern std::atomic<size_t> OMPallocator_device_mem_allocated;

inline size_t getOMPdeviceMemAllocated() { return OMPallocator_device_mem_allocated; }

template<typename T>
T* getOffloadDevicePtr(T* host_ptr)
{
  T* device_ptr;
  PRAGMA_OFFLOAD("omp target data use_device_ptr(host_ptr)")
  {
    device_ptr = host_ptr;
  }
  return device_ptr;
}

template<typename T, class HostAllocator = std::allocator<T>>
struct OMPallocator : public HostAllocator
{
  using value_type    = typename HostAllocator::value_type;
  using size_type     = typename HostAllocator::size_type;
  using pointer       = typename HostAllocator::pointer;
  using const_pointer = typename HostAllocator::const_pointer;

  OMPallocator() = default;
  OMPallocator(const OMPallocator&) : device_ptr(nullptr) {}
  OMPallocator& operator=(const OMPallocator&) { device_ptr = nullptr; }
  template<class U, class V>
  OMPallocator(const OMPallocator<U, V>&) : device_ptr(nullptr)
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
    OMPallocator_device_mem_allocated += n * sizeof(T);
    device_ptr = getOffloadDevicePtr(pt);
    return pt;
  }

  void deallocate(value_type* pt, std::size_t n)
  {
    PRAGMA_OFFLOAD("omp target exit data map(delete:pt[0:n])")
    OMPallocator_device_mem_allocated -= n * sizeof(T);
    HostAllocator::deallocate(pt, n);
  }

  T* getDevicePtr() { return device_ptr; }
  const T* getDevicePtr() const { return device_ptr; }

private:
  // pointee is on device.
  T* device_ptr = nullptr;
};

template<typename T, class HostAllocator>
struct allocator_traits<OMPallocator<T, HostAllocator>>
{
  static const bool is_host_accessible = true;
  static const bool is_dual_space = true;

  static void fill_n(T* ptr, size_t n, const T& value)
  {
    allocator_traits<HostAllocator>::fill_n(ptr, n, value);
    //PRAGMA_OFFLOAD("omp target update to(ptr[:n])")
  }
};

} // namespace qmcplusplus

#endif
