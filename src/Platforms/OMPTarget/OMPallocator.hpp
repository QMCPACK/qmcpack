//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file
 */
#ifndef QMCPLUSPLUS_OMPTARGET_ALLOCATOR_H
#define QMCPLUSPLUS_OMPTARGET_ALLOCATOR_H

#include <memory>
#include <type_traits>
#include <atomic>
#include "config.h"
#include "allocator_traits.hpp"
#ifdef ENABLE_CUDA
#include <cuda_runtime_api.h>
#endif
namespace qmcplusplus
{
extern std::atomic<size_t> OMPallocator_device_mem_allocated;

inline size_t getOMPdeviceMemAllocated() { return OMPallocator_device_mem_allocated; }

template<typename T>
T* getOffloadDevicePtr(T* host_ptr)
{
  T* device_ptr;
  PRAGMA_OFFLOAD("omp target data use_device_ptr(host_ptr)") { device_ptr = host_ptr; }
  return device_ptr;
}

/** OMPallocator is an allocator with fused device and dualspace allocator functionality.
 *  it is mostly c++03 style but is stateful with respect to the bond between the returned pt on the host
 *  and the device_ptr_.  While many containers may need a copy only one can own the memory
 *  it returns and it can only service one owner. i.e. only one object should call the allocate
 *  and deallocate methods.
 */
template<typename T, class HostAllocator = std::allocator<T>>
struct OMPallocator : public HostAllocator
{
  using value_type    = typename HostAllocator::value_type;
  using size_type     = typename HostAllocator::size_type;
  using pointer       = typename HostAllocator::pointer;
  using const_pointer = typename HostAllocator::const_pointer;

  OMPallocator() = default;
  /** Gives you a OMPallocator with no state.
   *  But OMPallocoator is stateful so this copy constructor is a lie.
   *  However until allocators are correct > c++11 this is retained since
   *  our < c++11 compliant containers may expect it.
   */
  OMPallocator(const OMPallocator&) : device_ptr_(nullptr) {}
  /** Used as a copy constructor when we need to actually copy the allocator.
   */
  OMPallocator(T* device_ptr, T* allocator_host_ptr) : device_ptr_(device_ptr), allocator_host_ptr_(allocator_host_ptr)
  {}

  // these semantics are surprising and "incorrect" considering OMPallocator is stateful.
  OMPallocator& operator=(const OMPallocator&) { device_ptr_ = nullptr; }

  // Use this if you actually expect need to assign.
  void actually_copy(const OMPallocator& from)
  {
    device_ptr_         = from.device_ptr_;
    allocator_host_ptr_ = from.allocator_host_ptr_;
  }

  template<class U, class V>
  OMPallocator(const OMPallocator<U, V>&) : device_ptr_(nullptr)
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
    device_ptr_         = getOffloadDevicePtr(pt);
    allocator_host_ptr_ = pt;
    return pt;
  }

  void deallocate(value_type* pt, std::size_t n)
  {
    PRAGMA_OFFLOAD("omp target exit data map(delete:pt[0:n])")
    OMPallocator_device_mem_allocated -= n * sizeof(T);
    HostAllocator::deallocate(pt, n);
  }

  T* getDevicePtr() { return device_ptr_; }
  const T* getDevicePtr() const { return device_ptr_; }

  T* getDevicePtr(T* host_ptr) { return device_ptr_ + (host_ptr - allocator_host_ptr_); }
  const T* getDevicePtr(T* host_ptr) const { return device_ptr_ + (host_ptr - allocator_host_ptr_); }

private:
  // pointee is on device.
  T* device_ptr_         = nullptr;
  T* allocator_host_ptr_ = nullptr;

public:
  T* get_device_ptr() const { return device_ptr_; }
  T* get_allocator_host_ptr() const { return allocator_host_ptr_; }
};

/** Specialization for OMPallocator which is a special DualAllocator with fused
 *  device and dualspace allocator functionality.
 */
template<typename T, class HostAllocator>
struct qmc_allocator_traits<OMPallocator<T, HostAllocator>>
{
  static const bool is_host_accessible = true;
  static const bool is_dual_space      = true;

  static void fill_n(T* ptr, size_t n, const T& value)
  {
    qmc_allocator_traits<HostAllocator>::fill_n(ptr, n, value);
    //PRAGMA_OFFLOAD("omp target update to(ptr[:n])")
  }

  static void updateTo(OMPallocator<T, HostAllocator>& alloc, T* host_ptr, size_t n)
  {
    T* device_ptr = alloc.getDevicePtr(host_ptr);
    #if defined(ENABLE_CUDA)
    cudaMemcpyAsync(device_ptr, host_ptr, sizeof(T) * n, cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();
    #endif
    // I want this to do what the above does, but obviously don't understand target udpate"
    PRAGMA_OFFLOAD("omp target update to(device_ptr[0:n])");
  }

  static void updateFrom(OMPallocator<T, HostAllocator>& alloc, T* host_ptr, size_t n)
  {
    T* device_ptr = alloc.getDevicePtr(host_ptr);
#if defined(ENABLE_CUDA)
    cudaMemcpyAsync(host_ptr, device_ptr, sizeof(T) * n, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
#endif
    // I want this to do what the above does, but obviously don't understand target udpate"
    PRAGMA_OFFLOAD("omp target update from(device_ptr[0:n])");
  }
};

} // namespace qmcplusplus
#endif
