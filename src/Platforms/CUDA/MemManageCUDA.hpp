//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file MemManageCUDA.hpp
 * this file provides three C++ memory allocators using CUDA specific memory allocation functions.
 *
 * CUDAManagedAllocator allocates CUDA unified memory
 * DeviceAllocator allocates CUDA device memory
 * CUDAHostAllocator allocates CUDA host pinned memory
 */
#ifndef QMCPLUSPLUS_MEMMANAGE_CUDA_H
#define QMCPLUSPLUS_MEMMANAGE_CUDA_H

#include <memory>
#include <cstdlib>
#include <stdexcept>
#include <atomic>
#include <limits>
#include <Common/MemManage.hpp>
#include "QueueCUDA.hpp"
#include "CUDAruntime.hpp"
#include "allocator_traits.hpp"
#include "CUDAfill.hpp"

namespace qmcplusplus
{

namespace compute
{

template<>
class MemManage<PlatformKind::CUDA>
{
public:
  static std::atomic<size_t> device_mem_allocated_;

  static size_t getDeviceMemAllocated() { return device_mem_allocated_; }

  static size_t getDeviceFreeMem()
  {
    size_t free, total;
    cudaErrorCheck(cudaMemGetInfo(&free, &total), "cudaMemGetInfo failed!");
    return free;
  }

  static void mallocDevice(void** ptr, size_t size) { cudaErrorCheck(cudaMalloc(ptr, size), "cudaMalloc failed!"); }

  static void freeDevice(void* ptr) { cudaErrorCheck(cudaFree(ptr), "cudaFree failed!"); }

  static void mallocHost(void** ptr, size_t size)
  {
    cudaErrorCheck(cudaMallocHost(ptr, size), "cudaMallocHost failed!");
  }

  static void freeHost(void* ptr) { cudaErrorCheck(cudaFreeHost(ptr), "cudaFreeHost failed!"); }

  static void registerHost(void* ptr, size_t size)
  {
    cudaErrorCheck(cudaHostRegister(ptr, size, cudaHostRegisterDefault), "cudaHostRegister failed!");
  }

  static void unregisterHost(void* ptr) { cudaErrorCheck(cudaHostUnregister(ptr), "cudaHostUnregister failed!"); }

  static void memcpy(void* dst, const void* src, size_t size)
  {
    cudaErrorCheck(cudaMemcpy(dst, src, size, cudaMemcpyDefault), "cudaMemcpy failed");
  }

  /// allocator for CUDA device memory
  template<typename T>
  using DeviceAllocator = DeviceAllocatorImpl<T, PlatformKind::CUDA>;

  /// allocator for CUDA host memory
  template<typename T>
  using HostAllocator = HostAllocatorImpl<T, PlatformKind::CUDA>;

  /// allocator for CUDA page locked memory
  template<typename T, class ULPHA = std::allocator<T>>
  using PageLockedAllocator = PageLockedAllocatorImpl<T, PlatformKind::CUDA, ULPHA>;
};

extern template class MemManage<PlatformKind::CUDA>;

} // namespace compute

/** allocator for CUDA unified memory
 * @tparam T data type
 */
template<typename T>
struct CUDAManagedAllocator
{
  using value_type    = T;
  using size_type     = size_t;
  using pointer       = T*;
  using const_pointer = const T*;

  CUDAManagedAllocator() = default;
  template<class U>
  CUDAManagedAllocator(const CUDAManagedAllocator<U>&)
  {}

  template<class U>
  struct rebind
  {
    using other = CUDAManagedAllocator<U>;
  };

  T* allocate(std::size_t n)
  {
    void* pt;
    cudaErrorCheck(cudaMallocManaged(&pt, n * sizeof(T)), "Allocation failed in CUDAManagedAllocator!");
    if ((size_t(pt)) & (QMC_SIMD_ALIGNMENT - 1))
      throw std::runtime_error("Unaligned memory allocated in CUDAManagedAllocator");
    return static_cast<T*>(pt);
  }
  void deallocate(T* p, std::size_t) { cudaErrorCheck(cudaFree(p), "Deallocation failed in CUDAManagedAllocator!"); }
};

template<class T>
using CUDAAllocator = compute::MemManage<PlatformKind::CUDA>::DeviceAllocator<T>;

template<typename T>
struct qmc_allocator_traits<qmcplusplus::CUDAAllocator<T>>
{
  static const bool is_host_accessible = false;
  static const bool is_dual_space      = false;
  static void fill_n(T* ptr, size_t n, const T& value) { qmcplusplus::CUDAfill_n(ptr, n, value); }
};

template<typename T>
using CUDAHostAllocator = compute::MemManage<PlatformKind::CUDA>::HostAllocator<T>;

template<typename T, class ULPHA = std::allocator<T>>
using CUDALockedPageAllocator = compute::MemManage<PlatformKind::CUDA>::PageLockedAllocator<T, ULPHA>;

} // namespace qmcplusplus

#endif
