//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file MemManageSYCL.hpp
 * this file provides three C++ memory allocators using SYCL specific memory allocation functions.
 *
 * SYCLManagedAllocator allocates SYCL shared memory
 * DeviceAllocator allocates SYCL device memory
 * SYCLHostAllocator allocates SYCL host memory
 * They are based on CUDA*Allocator implementation
 */
#ifndef QMCPLUSPLUS_MEMMANAGE_SYCL_H
#define QMCPLUSPLUS_MEMMANAGE_SYCL_H

#include <memory>
#include <cstdlib>
#include <stdexcept>
#include <atomic>
#include <limits>
#include "config.h"
#include <Common/MemManage.hpp>
#include "QueueSYCL.hpp"
#include "allocator_traits.hpp"
#include "SYCLruntime.hpp"

namespace qmcplusplus
{
extern std::atomic<size_t> SYCLallocator_device_mem_allocated;

inline size_t getSYCLdeviceMemAllocated() { return SYCLallocator_device_mem_allocated; }

namespace compute
{

template<>
class MemManage<PlatformKind::SYCL>
{
public:
  static void mallocDevice(void** ptr, size_t size)
  {
    *ptr = sycl::malloc_device(size, getSYCLDefaultDeviceDefaultQueue());
  }

  static void freeDevice(void* ptr) { sycl::free(ptr, getSYCLDefaultDeviceDefaultQueue()); }

  static void registerHost(void* ptr, size_t size)
  {
    sycl::ext::oneapi::experimental::prepare_for_device_copy(ptr, size, getSYCLDefaultDeviceDefaultQueue());
  }

  static void unregisterHost(void* ptr)
  {
    sycl::ext::oneapi::experimental::release_from_device_copy(ptr, getSYCLDefaultDeviceDefaultQueue());
  }

  /** allocator for SYCL device memory
   * @tparm T data type
   * @tparm ALIGN alignment in bytes
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
      MemManage<PlatformKind::SYCL>::mallocDevice(&pt, n * sizeof(T));
      SYCLallocator_device_mem_allocated += n * sizeof(T);
      return static_cast<T*>(pt);
    }

    void deallocate(T* p, std::size_t n)
    {
      MemManage<PlatformKind::SYCL>::freeDevice(p);
      SYCLallocator_device_mem_allocated -= n * sizeof(T);
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

    void copyToDevice(T* device_ptr, T* host_ptr, size_t n)
    {
      getSYCLDefaultDeviceDefaultQueue().memcpy(device_ptr, host_ptr, n * sizeof(T)).wait();
    }

    void copyFromDevice(T* host_ptr, T* device_ptr, size_t n)
    {
      getSYCLDefaultDeviceDefaultQueue().memcpy(host_ptr, device_ptr, n * sizeof(T)).wait();
    }

    void copyDeviceToDevice(T* to_ptr, size_t n, T* from_ptr)
    {
      getSYCLDefaultDeviceDefaultQueue().memcpy(to_ptr, from_ptr, n * sizeof(T)).wait();
    }
  };

  /** allocator for SYCL host pinned memory
 * @tparm T data type
 * @tparm ALIGN alignment in bytes
 */
  template<typename T, size_t ALIGN = QMC_SIMD_ALIGNMENT>
  struct PinnedAllocator
  {
    using value_type    = T;
    using size_type     = size_t;
    using pointer       = T*;
    using const_pointer = const T*;

    static constexpr size_t alignment = ALIGN;

    PinnedAllocator() = default;
    template<class U>
    PinnedAllocator(const PinnedAllocator<U>&)
    {}

    template<class U>
    struct rebind
    {
      typedef PinnedAllocator<U, ALIGN> other;
    };

    T* allocate(std::size_t n) { return sycl::aligned_alloc_host<T>(ALIGN, n, getSYCLDefaultDeviceDefaultQueue()); }
    void deallocate(T* p, std::size_t) { sycl::free(p, getSYCLDefaultDeviceDefaultQueue()); }
  };
};
extern template class MemManage<PlatformKind::SYCL>;
}; // namespace compute


/** allocator for SYCL shared memory
 * @tparm T data type
 * @tparm ALIGN alignment in bytes
 */
template<typename T>
struct SYCLSharedAllocator
{
  typedef T value_type;
  typedef size_t size_type;
  typedef T* pointer;
  typedef const T* const_pointer;

  SYCLSharedAllocator() = default;
  template<class U>
  SYCLSharedAllocator(const SYCLSharedAllocator<U>&)
  {}

  template<class U>
  struct rebind
  {
    typedef SYCLSharedAllocator<U> other;
  };

  T* allocate(std::size_t n)
  {
    T* pt = sycl::malloc_shared<T>(n, getSYCLDefaultDeviceDefaultQueue());
    return pt;
  }

  void deallocate(T* p, std::size_t) { sycl::free(p, getSYCLDefaultDeviceDefaultQueue()); }
};

template<class T>
using SYCLAllocator = compute::MemManage<PlatformKind::SYCL>::DeviceAllocator<T>;

template<typename T>
struct qmc_allocator_traits<qmcplusplus::SYCLAllocator<T>>
{
  static const bool is_host_accessible = false;
  static const bool is_dual_space      = false;
  static void fill_n(T* ptr, size_t n, const T& value)
  {
    //THINK
    //qmcplusplus::SYCLfill_n(ptr, n, value);
  }
  static void updateTo(SYCLAllocator<T>& alloc, T* host_ptr, size_t n)
  {
    T* device_ptr = alloc.getDevicePtr(host_ptr);
    alloc.copyToDevice(device_ptr, host_ptr, n);
  }

  static void updateFrom(SYCLAllocator<T>& alloc, T* host_ptr, size_t n)
  {
    T* device_ptr = alloc.getDevicePtr(host_ptr);
    alloc.copyFromDevice(host_ptr, device_ptr, n);
  }
};

template<typename T, size_t ALIGN = QMC_SIMD_ALIGNMENT>
using SYCLHostAllocator = compute::MemManage<PlatformKind::SYCL>::PinnedAllocator<T, ALIGN>;
} // namespace qmcplusplus

#endif
