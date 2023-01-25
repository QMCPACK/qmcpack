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
/** @file SYCLallocator.hpp
 * this file provides three C++ memory allocators using SYCL specific memory allocation functions.
 *
 * SYCLManagedAllocator allocates SYCL shared memory
 * SYCLAllocator allocates SYCL device memory
 * SYCLHostAllocator allocates SYCL host memory
 * They are based on CUDA*Allocator implementation
 */
#ifndef QMCPLUSPLUS_SYCL_ALLOCATOR_H
#define QMCPLUSPLUS_SYCL_ALLOCATOR_H

#include <memory>
#include <cstdlib>
#include <stdexcept>
#include <atomic>
#include <limits>
#include <sycl/sycl.hpp>
#include "config.h"
#include "allocator_traits.hpp"
#include "SYCLruntime.hpp"

namespace qmcplusplus
{
extern std::atomic<size_t> SYCLallocator_device_mem_allocated;

inline size_t getSYCLdeviceMemAllocated() { return SYCLallocator_device_mem_allocated; }

/** allocator for SYCL shared memory
 * @tparm T data type
 * @tparm ALIGN alignment in bytes
 */
template<typename T, size_t ALIGN = QMC_SIMD_ALIGNMENT>
struct SYCLSharedAllocator
{
  typedef T value_type;
  typedef size_t size_type;
  typedef T* pointer;
  typedef const T* const_pointer;

  static constexpr size_t alignment = ALIGN;

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
    T* pt = sycl::aligned_alloc_shared<T>(ALIGN, n, getSYCLDefaultDeviceDefaultQueue());
    return pt;
  }
  void deallocate(T* p, std::size_t) { sycl::free(p, getSYCLDefaultDeviceDefaultQueue()); }
};

template<class T1, class T2>
bool operator==(const SYCLSharedAllocator<T1>&, const SYCLSharedAllocator<T2>&)
{
  return true;
}
template<class T1, class T2>
bool operator!=(const SYCLSharedAllocator<T1>&, const SYCLSharedAllocator<T2>&)
{
  return false;
}

/** allocator for SYCL device memory
 * @tparm T data type
 * @tparm ALIGN alignment in bytes
 *
 * using this with something other than Ohmms containers?
 *  -- use caution, write unit tests! --
 * It's not tested beyond use in some unit tests using std::vector with constant size.
 * SYCLAllocator appears to meet all the nonoptional requirements of a c++ Allocator.
 *
 * Some of the default implementations in std::allocator_traits
 * of optional Allocator requirements may cause runtime or compilation failures.
 * They assume there is only one memory space and that the host has access to it.
 */
template<typename T, size_t ALIGN = 64>
class SYCLAllocator
{
public:
  typedef T value_type;
  typedef size_t size_type;
  typedef T* pointer;
  typedef const T* const_pointer;

  static constexpr size_t alignment = ALIGN;

  SYCLAllocator() = default;
  template<class U>
  SYCLAllocator(const SYCLAllocator<U>&)
  {}

  template<class U>
  struct rebind
  {
    typedef SYCLAllocator<U> other;
  };

  T* allocate(std::size_t n)
  {
    T* pt = sycl::aligned_alloc_device<T>(ALIGN, n, getSYCLDefaultDeviceDefaultQueue());
    SYCLallocator_device_mem_allocated += n * sizeof(T);
    return pt;
  }

  void deallocate(T* p, std::size_t n)
  {
    sycl::free(p, getSYCLDefaultDeviceDefaultQueue());
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

template<class T1, class T2>
bool operator==(const SYCLAllocator<T1>&, const SYCLAllocator<T2>&)
{
  return true;
}
template<class T1, class T2>
bool operator!=(const SYCLAllocator<T1>&, const SYCLAllocator<T2>&)
{
  return false;
}

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

/** allocator for SYCL host pinned memory
 * @tparm T data type
 * @tparm ALIGN alignment in bytes
 */
template<typename T, size_t ALIGN = QMC_SIMD_ALIGNMENT>
struct SYCLHostAllocator
{
  typedef T value_type;
  typedef size_t size_type;
  typedef T* pointer;
  typedef const T* const_pointer;

  static constexpr size_t alignment = ALIGN;

  SYCLHostAllocator() = default;
  template<class U>
  SYCLHostAllocator(const SYCLHostAllocator<U>&)
  {}

  template<class U>
  struct rebind
  {
    typedef SYCLHostAllocator<U> other;
  };

  T* allocate(std::size_t n) { return sycl::aligned_alloc_host<T>(ALIGN, n, getSYCLDefaultDeviceDefaultQueue()); }
  void deallocate(T* p, std::size_t) { sycl::free(p, getSYCLDefaultDeviceDefaultQueue()); }
};

template<class T1, class T2>
bool operator==(const SYCLHostAllocator<T1>&, const SYCLHostAllocator<T2>&)
{
  return true;
}

template<class T1, class T2>
bool operator!=(const SYCLHostAllocator<T1>&, const SYCLHostAllocator<T2>&)
{
  return false;
}

} // namespace qmcplusplus

#endif
