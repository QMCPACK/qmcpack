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
  static std::atomic<size_t> device_mem_allocated_;

public:
  static size_t getDeviceMemAllocated() { return device_mem_allocated_; }

  static size_t getDeviceFreeMem()
  {
    size_t free, total;
    cudaErrorCheck(cudaMemGetInfo(&free, &total), "cudaMemGetInfo failed!");
    return free;
  }
  static void mallocDevice(void** ptr, size_t size) { cudaErrorCheck(cudaMalloc(ptr, size), "cudaMalloc failed!"); }

  static void freeDevice(void* ptr) { cudaErrorCheck(cudaFree(ptr), "cudaFree failed!"); }

  static void registerHost(void* ptr, size_t size)
  {
    cudaErrorCheck(cudaHostRegister(ptr, size, cudaHostRegisterDefault), "cudaHostRegister failed!");
  }

  static void unregisterHost(void* ptr) { cudaErrorCheck(cudaHostUnregister(ptr), "cudaHostUnregister failed!"); }

  static void memcpy(void* dst, const void* src, size_t size)
  {
    cudaErrorCheck(cudaMemcpy(dst, src, size, cudaMemcpyDefault), "cudaMemcpy failed");
  }
  /** allocator for CUDA device memory
   * @tparam T data type
   *
   * using this with something other than Ohmms containers?
   *  -- use caution, write unit tests! --
   * It's not tested beyond use in some unit tests using std::vector with constant size.
   * DeviceAllocator appears to meet all the nonoptional requirements of a c++ Allocator.
   *
   * Some of the default implementations in std::allocator_traits
   * of optional Allocator requirements may cause runtime or compilation failures.
   * They assume there is only one memory space and that the host has access to it.
   */
  template<typename T>
  class DeviceAllocator
  {
  public:
    using value_type    = T;
    using size_type     = size_t;
    using pointer       = T*;
    using const_pointer = const T*;

    DeviceAllocator() = default;
    template<class U>
    DeviceAllocator(const DeviceAllocator<U>&)
    {}

    template<class U>
    struct rebind
    {
      using other = DeviceAllocator<U>;
    };

    T* allocate(std::size_t n)
    {
      void* pt;
      MemManage<PlatformKind::CUDA>::mallocDevice(&pt, n * sizeof(T));
      device_mem_allocated_ += n * sizeof(T);
      return static_cast<T*>(pt);
    }

    void deallocate(T* p, std::size_t n)
    {
      MemManage<PlatformKind::CUDA>::freeDevice(p);
      device_mem_allocated_ -= n * sizeof(T);
    }

    /** Provide a construct for std::allocator_traits::contruct to call.
     *  Don't do anything on construct, pointer p is on the device!
     *
     *  For example std::vector calls this to default initialize each element. You'll segfault
     *  if std::allocator_traits::construct tries doing that at p.
     *
     *  The standard is a bit confusing on this point. Implementing this is an optional requirement
     *  of Allocator from C++11 on, its not slated to be removed.
     *
     *  Its deprecated for the std::allocator in c++17 and will be removed in c++20.  But we are not implementing
     *  std::allocator.
     *
     *  STL containers only use Allocators through allocator_traits and std::allocator_traits handles the case
     *  where no construct method is present in the Allocator.
     *  But std::allocator_traits will call the Allocators construct method if present.
     */
    template<class U, class... Args>
    static void construct(U* p, Args&&... args)
    {}

    /** Give std::allocator_traits something to call.
     *  The default if this isn't present is to call p->~T() which
     *  we can't do on device memory.
     */
    template<class U>
    static void destroy(U* p)
    {}

    static void memcpy(T* dst, const T* src, size_t size)
    {
      MemManage<PlatformKind::CUDA>::memcpy(dst, src, sizeof(T) * size);
    }
  };

  /** allocator locks memory pages allocated by ULPHA
   * @tparam T data type
   * @tparam ULPHA host memory allocator using unlocked page
   *
   * ULPHA cannot be CUDAHostAllocator
   */
  template<typename T, class ULPHA = std::allocator<T>>
  struct PinnedAllocator : public ULPHA
  {
    using value_type    = typename ULPHA::value_type;
    using size_type     = typename ULPHA::size_type;
    using pointer       = typename ULPHA::pointer;
    using const_pointer = typename ULPHA::const_pointer;

    PinnedAllocator() = default;
    template<class U, class V>
    PinnedAllocator(const PinnedAllocator<U, V>&)
    {}

    template<class U>
    struct rebind
    {
      using other = PinnedAllocator<U, typename std::allocator_traits<ULPHA>::template rebind_alloc<U>>;
    };

    value_type* allocate(std::size_t n)
    {
      static_assert(std::is_same<T, value_type>::value, "PinnedAllocator and ULPHA data types must agree!");
      value_type* pt = ULPHA::allocate(n);
      MemManage<PlatformKind::CUDA>::registerHost(pt, n * sizeof(T));
      return pt;
    }

    void deallocate(value_type* pt, std::size_t n)
    {
      MemManage<PlatformKind::CUDA>::unregisterHost(pt);
      ULPHA::deallocate(pt, n);
    }
  };
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
using CUDAHostAllocator = compute::MemManage<PlatformKind::CUDA>::PinnedAllocator<T>;

template<typename T, class ULPHA = std::allocator<T>>
using CUDALockedPageAllocator = compute::MemManage<PlatformKind::CUDA>::PinnedAllocator<T, ULPHA>;

} // namespace qmcplusplus

#endif
