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
namespace compute
{

template<>
class MemManage<PlatformKind::SYCL>
{
public:
  static std::atomic<size_t> device_mem_allocated_;

  static size_t getDeviceMemAllocated() { return device_mem_allocated_; }

  static size_t getDeviceFreeMem()
  {
    auto device = getSYCLDefaultDeviceDefaultQueue().get_device();
    if (device.has(sycl::aspect::ext_intel_free_memory))
      return getSYCLDefaultDeviceDefaultQueue().get_device().get_info<sycl::ext::intel::info::device::free_memory>();
    else
      return 0;
  }

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

  static void memcpy(void* dst, const void* src, size_t size)
  {
    getSYCLDefaultDeviceDefaultQueue().memcpy(dst, src, size).wait();
  }

  /// allocator for SYCL device memory
  template<typename T>
  using DeviceAllocator = DeviceAllocatorImpl<T, PlatformKind::SYCL>;

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
};

template<typename T, size_t ALIGN = QMC_SIMD_ALIGNMENT>
using SYCLHostAllocator = compute::MemManage<PlatformKind::SYCL>::PinnedAllocator<T, ALIGN>;
} // namespace qmcplusplus

#endif
