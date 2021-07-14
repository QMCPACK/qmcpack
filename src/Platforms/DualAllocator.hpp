//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// Refactored from: OMPallocator.hpp
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file
 */
#ifndef QMCPLUSPLUS_DUAL_ALLOCATOR_H
#define QMCPLUSPLUS_DUAL_ALLOCATOR_H

#include <memory>
#include <type_traits>
#include <atomic>
#include <exception>
#include "config.h"
#include "allocator_traits.hpp"
#include "PinnedAllocator.h"
#if defined(ENABLE_CUDA)
#include "CUDA/CUDAallocator.hpp"
#endif

namespace qmcplusplus
{
extern std::atomic<size_t> dual_device_mem_allocated;

inline size_t getDualDeviceMemAllocated() { return dual_device_mem_allocated; }

/** Generalizes the DualMemorySpace allocator
 *  This provides a limited alternative to OMPallocator for testing/benchmarking
 *  without dependence of OMPTarget/ offload.
 *  It does not provide an alternative to OMPtarget transfer semantics so many production
 *  objects will not be functional if it is used as the allocator for the data objects they depend
 *  on.
 *  If you use DualAllocator at this time you need to handle data transfer yourself.
 *
 *  \todo the OMPTarget allocation can be a "device" allocator comparable to CUDAAllocator
 *  Then OMPallocator can be replaced by a DualAllocator<T, OffloadAllocator<T>, PinnedAllocator<T>>
 */
template<typename T, class DeviceAllocator, class HostAllocator = std::allocator<T>>
struct DualAllocator : public HostAllocator
{
  using Value        = typename HostAllocator::value_type;
  using Size         = typename HostAllocator::size_type;
  using Pointer      = typename HostAllocator::pointer;
  using ConstPointer = typename HostAllocator::const_pointer;

  DualAllocator() : device_ptr_(nullptr) {};
  DualAllocator(const DualAllocator&) : device_ptr_(nullptr) {}
  DualAllocator& operator=(const DualAllocator&) { device_ptr_ = nullptr; }
  template<class U, class V>
  DualAllocator(const DualAllocator<U, V>&) : device_ptr_(nullptr)
  {}

  template<class U, class V>
  struct rebind
  {
    using other = DualAllocator<U, V>;
  };

  Value* allocate(std::size_t n)
  {
    static_assert(std::is_same<T, Value>::value, "DualAllocator and HostAllocator data types must agree!");
    if (device_ptr_ != nullptr)
      throw std::runtime_error("DualAllocator does not support device reallocation");
    T* host_ptr   = std::allocator_traits<HostAllocator>::allocate(allocator_, n);
    device_ptr_ = std::allocator_traits<DeviceAllocator>::allocate(device_allocator_, n);
    dual_device_mem_allocated += n * sizeof(T);
    return host_ptr;
  }

  void deallocate(Value* pt, std::size_t n)
  {
    dual_device_mem_allocated -= n * sizeof(T);
    std::allocator_traits<DeviceAllocator>::deallocate(device_allocator_, device_ptr_, n);
    std::allocator_traits<HostAllocator>::deallocate(allocator_, pt, n);
    device_ptr_ = nullptr;
  }

  void attachReference(DualAllocator& from, std::ptrdiff_t ptr_offset)
  {
    device_ptr_ = from.getDevicePtr() + ptr_offset;
  }

  T* getDevicePtr() { return device_ptr_; }
  const T* getDevicePtr() const { return device_ptr_; }

  DeviceAllocator& get_device_allocator() { return device_allocator_; }
  const DeviceAllocator& get_device_allocator() const { return device_allocator_; }

private:
  HostAllocator allocator_;
  DeviceAllocator device_allocator_;
  T* device_ptr_;
  // host pointer to device pointer map, needed for deallocation
};

template<typename T, class DeviceAllocator, class HostAllocator>
struct qmc_allocator_traits<DualAllocator<T, DeviceAllocator, HostAllocator>>
{
  using DualAlloc = DualAllocator<T, DeviceAllocator, HostAllocator>;
  static const bool is_host_accessible = true;
  static const bool is_dual_space      = true;

  static void fill_n(T* ptr, size_t n, const T& value) { qmc_allocator_traits<HostAllocator>::fill_n(ptr, n, value); }

  static void attachReference(DualAlloc& from, DualAlloc& to, std::ptrdiff_t ptr_offset)
  {
    to.attachReference(from, ptr_offset);
  }

  static void updateTo(DualAlloc& alloc, T* host_ptr, size_t n)
  {
    alloc.get_device_allocator().copyToDevice(alloc.getDevicePtr(), host_ptr, n);
  }

  static void updateFrom(DualAlloc& alloc, T* host_ptr, size_t n)
  {
    alloc.get_device_allocator().copyFromDevice(host_ptr, alloc.getDevicePtr(), n);
  }

  static void deviceSideCopyN(DualAlloc& alloc, size_t to, size_t n, size_t from) {
    T* device_ptr = alloc.getDevicePtr();
    T* to_ptr = device_ptr + to;
    T* from_ptr = device_ptr + from;
    alloc.get_device_allocator().copyDeviceToDevice(to_ptr, n, from_ptr);
  }
};

} // namespace qmcplusplus

#endif
