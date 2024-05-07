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
#if defined(ENABLE_OFFLOAD)
#include <omp.h>
#endif

#if defined(QMC_OFFLOAD_MEM_ASSOCIATED)
#include <CUDA/CUDAruntime.hpp>
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
 *
 *  Note: in the style of openmp portability this class always thinks its dual space even when its not,
 *  this happens through the magic of openmp ignoring target pragmas when target isn't enabled. i.e.
 *  -fopenmp-targets=... isn't passed at compile time.
 *  This makes the code "simpler" and more "portable" since its the same code you would write for
 *  openmp CPU implementation *exploding head* and that is the same implementation + pragmas
 *  as the serial implementation. This definitely isn't true for all QMCPACK code using offload
 *  but it is true for OMPAllocator so we do test it that way.
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
  template<class U, class V>
  OMPallocator(const OMPallocator<U, V>&) : device_ptr_(nullptr)
  {}

  template<class U, class V>
  struct rebind
  {
    using other = OMPallocator<U, V>;
  };

  value_type* allocate(std::size_t n)
  {
    static_assert(std::is_same<T, value_type>::value, "OMPallocator and HostAllocator data types must agree!");
    value_type* pt = HostAllocator::allocate(n);
#if defined(QMC_OFFLOAD_MEM_ASSOCIATED)
    cudaErrorCheck(cudaMalloc(&device_ptr_, n * sizeof(T)), "cudaMalloc failed in OMPallocator!");
    const int status = omp_target_associate_ptr(pt, device_ptr_, n * sizeof(T), 0, omp_get_default_device());
    if (status != 0)
      throw std::runtime_error("omp_target_associate_ptr failed in OMPallocator!");
#else
    PRAGMA_OFFLOAD("omp target enter data map(alloc:pt[0:n])")
    device_ptr_ = getOffloadDevicePtr(pt);
#endif
    OMPallocator_device_mem_allocated += n * sizeof(T);
    return pt;
  }

  void deallocate(value_type* pt, std::size_t n)
  {
    OMPallocator_device_mem_allocated -= n * sizeof(T);
#if defined(QMC_OFFLOAD_MEM_ASSOCIATED)
    T* device_ptr_from_omp = getOffloadDevicePtr(pt);
    const int status       = omp_target_disassociate_ptr(pt, omp_get_default_device());
    if (status != 0)
      throw std::runtime_error("omp_target_disassociate_ptr failed in OMPallocator!");
    cudaErrorCheck(cudaFree(device_ptr_from_omp), "cudaFree failed in OMPallocator!");
#else
    PRAGMA_OFFLOAD("omp target exit data map(delete:pt[0:n])")
#endif
    HostAllocator::deallocate(pt, n);
  }

  void attachReference(const OMPallocator& from, std::ptrdiff_t ptr_offset)
  {
    device_ptr_ = const_cast<typename OMPallocator::pointer>(from.get_device_ptr()) + ptr_offset;
  }

  T* get_device_ptr() { return device_ptr_; }
  const T* get_device_ptr() const { return device_ptr_; }

private:
  // pointee is on device.
  T* device_ptr_ = nullptr;
};

/** Specialization for OMPallocator which is a special DualAllocator with fused
 *  device and dualspace allocator functionality.
 */
template<typename T, class HostAllocator>
struct qmc_allocator_traits<OMPallocator<T, HostAllocator>>
{
  static constexpr bool is_host_accessible = true;
  static constexpr bool is_dual_space      = true;

  static void fill_n(T* ptr, size_t n, const T& value)
  {
    qmc_allocator_traits<HostAllocator>::fill_n(ptr, n, value);
    //PRAGMA_OFFLOAD("omp target update to(ptr[:n])")
  }

  static void attachReference(const OMPallocator<T, HostAllocator>& from,
                              OMPallocator<T, HostAllocator>& to,
                              std::ptrdiff_t ptr_offset)
  {
    to.attachReference(from, ptr_offset);
  }

  static void updateTo(OMPallocator<T, HostAllocator>& alloc, T* host_ptr, size_t n, size_t offset = 0)
  {
    PRAGMA_OFFLOAD("omp target update to(host_ptr[offset:n])");
  }

  static void updateFrom(OMPallocator<T, HostAllocator>& alloc, T* host_ptr, size_t n, size_t offset = 0)
  {
    PRAGMA_OFFLOAD("omp target update from(host_ptr[offset:n])");
  }

  // Not very optimized device side copy.  Only used for testing.
  static void deviceSideCopyN(OMPallocator<T, HostAllocator>& alloc, size_t to, size_t n, size_t from)
  {
    auto* dev_ptr = alloc.get_device_ptr();
    PRAGMA_OFFLOAD("omp target teams distribute parallel for is_device_ptr(dev_ptr)")
    for (int i = 0; i < n; i++)
      dev_ptr[to + i] = dev_ptr[from + i];
  }
};

#if defined(ENABLE_OFFLOAD)
/** allocator for OMPTarget device memory
 * @tparam T data type
 *
 * using this with something other than Ohmms containers?
 *  -- use caution, write unit tests! --
 * It's not tested beyond use in some unit tests using std::vector with constant size.
 * OMPTargetAllocator appears to meet all the nonoptional requirements of a c++ Allocator.
 *
 * Some of the default implementations in std::allocator_traits
 * of optional Allocator requirements may cause runtime or compilation failures.
 * They assume there is only one memory space and that the host has access to it.
 */
template<typename T>
class OMPTargetAllocator
{
public:
  using value_type    = T;
  using size_type     = size_t;
  using pointer       = T*;
  using const_pointer = const T*;

  OMPTargetAllocator() = default;
  template<class U>
  OMPTargetAllocator(const OMPTargetAllocator<U>&)
  {}

  template<class U>
  struct rebind
  {
    using other = OMPTargetAllocator<U>;
  };

  T* allocate(std::size_t n)
  {
    void* pt = omp_target_alloc(n * sizeof(T), omp_get_default_device());
    if(!pt)
      throw std::runtime_error("Allocation failed in OMPTargetAllocator!");
    OMPallocator_device_mem_allocated += n * sizeof(T);
    return static_cast<T*>(pt);
  }

  void deallocate(T* p, std::size_t n)
  {
    omp_target_free(p, omp_get_default_device());
    OMPallocator_device_mem_allocated -= n * sizeof(T);
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
    if(omp_target_memcpy(device_ptr, host_ptr, n, 0, 0, omp_get_default_device(), omp_get_initial_device()))
      throw std::runtime_error("omp_target_memcpy failed in copyToDevice");
  }

  void copyFromDevice(T* host_ptr, T* device_ptr, size_t n)
  {
    if(omp_target_memcpy(host_ptr, device_ptr, n, 0, 0, omp_get_initial_device(), omp_get_default_device()))
      throw std::runtime_error("omp_target_memcpy failed in copyToDevice");
  }

  void copyDeviceToDevice(T* to_ptr, size_t n, T* from_ptr)
  {
    if(omp_target_memcpy(to_ptr, from_ptr, n, 0, 0, omp_get_default_device(), omp_get_default_device()))
      throw std::runtime_error("omp_target_memcpy failed in copyToDevice");
  }
};

template<class T1, class T2>
bool operator==(const OMPTargetAllocator<T1>&, const OMPTargetAllocator<T2>&)
{
  return true;
}
template<class T1, class T2>
bool operator!=(const OMPTargetAllocator<T1>&, const OMPTargetAllocator<T2>&)
{
  return false;
}

template<typename T>
struct qmc_allocator_traits<qmcplusplus::OMPTargetAllocator<T>>
{
  static const bool is_host_accessible = false;
  static const bool is_dual_space      = false;
  static void fill_n(T* ptr, size_t n, const T& value) { }
};
#endif

} // namespace qmcplusplus
#endif
