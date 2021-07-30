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
#include "Synchro.hpp"
#include "CUDA/CUDAallocator.hpp"

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
  using Synchro_t      = typename DeviceAllocator::Synchro_t;
  
  DualAllocator() : host_ptr_(nullptr), device_ptr_(nullptr), synchro_(default_synchro) {};
  DualAllocator(const DualAllocator&) : host_ptr_(nullptr), device_ptr_(nullptr) {}
  DualAllocator& operator=(const DualAllocator&)
  {
    host_ptr_   = nullptr;
    device_ptr_ = nullptr;
  }
  template<class U, class V>
  DualAllocator(const DualAllocator<U, V>&) : host_ptr_(nullptr), device_ptr_(nullptr)
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
    host_ptr_   = std::allocator_traits<HostAllocator>::allocate(allocator_, n);
    device_ptr_ = std::allocator_traits<DeviceAllocator>::allocate(device_allocator_, n);
    dual_device_mem_allocated += n * sizeof(T);
    return host_ptr_;
  }

  void deallocate(Value* pt, std::size_t n)
  {
    assert(pt = host_ptr_);
    dual_device_mem_allocated -= n * sizeof(T);
    std::allocator_traits<DeviceAllocator>::deallocate(device_allocator_, device_ptr_, n);
    std::allocator_traits<HostAllocator>::deallocate(allocator_, pt, n);
    device_ptr_ = nullptr;
    host_ptr_   = nullptr;
  }

  void attachReference(const DualAllocator& from, const T* from_data, T* ref)
  {
    std::ptrdiff_t ptr_offset = ref - from_data;
    host_ptr_                 = ref;
    device_ptr_               = const_cast<typename DualAllocator::pointer>(from.get_device_ptr()) + ptr_offset;
  }

  T* get_device_ptr() { return device_ptr_; }
  const T* get_device_ptr() const { return device_ptr_; }

  DeviceAllocator& get_device_allocator() { return device_allocator_; }
  const DeviceAllocator& get_device_allocator() const { return device_allocator_; }

  T* get_host_ptr() { return host_ptr_; }

  template<class SYNCHRO>
  void setSync(SYNCHRO& synchro) { synchro_ = synchro; }  
  Synchro& getSync() { return synchro_; }
  void sync() { synchro_.get().sync(); }
private:
  HostAllocator allocator_;
  DeviceAllocator device_allocator_;
  T* host_ptr_;
  T* device_ptr_;
  // would rather std::optional I think.
  Synchro_t default_synchro;
  std::reference_wrapper<Synchro> synchro_;
};

template<typename T, class DeviceAllocator, class HostAllocator>
struct qmc_allocator_traits<DualAllocator<T, DeviceAllocator, HostAllocator>>
{
  using DualAlloc                      = DualAllocator<T, DeviceAllocator, HostAllocator>;
  static const bool is_host_accessible = true;
  static const bool is_dual_space      = true;

  static void fill_n(T* ptr, size_t n, const T& value) { qmc_allocator_traits<HostAllocator>::fill_n(ptr, n, value); }

  static void attachReference(const DualAlloc& from, DualAlloc& to, const T* from_data, T* ref)
  {
    to.attachReference(from, from_data, ref);
  }

  static void setSync(DualAlloc& alloc, Synchro& synchro) { alloc.setSync(synchro); }

  static void sync(DualAlloc& alloc, Synchro& synchro) { alloc.sync(); }

  /** update to the device, assumes you are copying starting with the implicit host_ptr.
   *
   *  These follow the openmp target semantics where you only provide the host
   *  side of a host_ptr device_ptr pair but the verb relates to what happens on the device.
   *
   *  This is primarily for testing to reduce ifdef code and single "flavor" testing
   *
   *  In a few places in the code where transfer crept up into otherwise "Device" abstracted code
   *  it is necessary.
   *
   *  This a generic API and unlikely to be the best way to handle performance critical transfers,
   *  but if you have to use it or ifdef at a level above a xxxCUDA.cu or xxxOMPTarget.hpp file
   *  thats an issue and this is the best answer if you can't factor it down into CUDA or OMPTarget specific
   *  code.
   */
  static void updateTo(DualAlloc& alloc, T* host_ptr, size_t n)
  {
    assert(host_ptr = alloc.get_host_ptr());
    alloc.get_device_allocator().copyToDevice(alloc.get_device_ptr(), host_ptr, n);
  }

  /** update from the device, assumes you are copying starting with the device_ptr to the implicit host_ptr.
   */
  static void updateFrom(DualAlloc& alloc, T* host_ptr, size_t n)
  {
    assert(host_ptr = alloc.get_host_ptr());
    alloc.get_device_allocator().copyFromDevice(host_ptr, alloc.get_device_ptr(), n);
  }

  static void updateToAsync(DualAlloc& alloc, T* host_ptr, size_t n)
  {
    assert(host_ptr = alloc.get_host_ptr());
    alloc.get_device_allocator().copyToDeviceAsync(alloc.get_device_ptr(), host_ptr, n, alloc.getSync());
  }

  /** update from the device, assumes you are copying starting with the device_ptr to the implicit host_ptr.
   */
  static void updateFromAsync(DualAlloc& alloc, T* host_ptr, size_t n)
  {
    assert(host_ptr = alloc.get_host_ptr());
    alloc.get_device_allocator().copyFromDeviceAsync(host_ptr, alloc.get_device_ptr(), n, alloc.getSync());
  }
  
  static void deviceSideCopyN(DualAlloc& alloc, size_t to, size_t n, size_t from)
  {
    T* device_ptr = alloc.get_device_ptr();
    T* to_ptr     = device_ptr + to;
    T* from_ptr   = device_ptr + from;
    alloc.get_device_allocator().copyDeviceToDevice(to_ptr, n, from_ptr);
  }
};

} // namespace qmcplusplus

#endif
