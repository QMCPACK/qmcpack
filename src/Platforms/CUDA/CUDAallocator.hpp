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
/** @file CUDAallocator.hpp
 * this file provides three C++ memory allocators using CUDA specific memory allocation functions.
 *
 * CUDAManagedAllocator allocates CUDA unified memory
 * CUDAAllocator allocates CUDA device memory
 * CUDAHostAllocator allocates CUDA host pinned memory
 */
#ifndef QMCPLUSPLUS_CUDA_ALLOCATOR_H
#define QMCPLUSPLUS_CUDA_ALLOCATOR_H

#include <memory>
#include <cstdlib>
#include <stdexcept>
#include <atomic>
#include <limits>
#include "CUDAruntime.hpp"
#include "allocator_traits.hpp"
#include "CUDAfill.hpp"

namespace qmcplusplus
{
extern std::atomic<size_t> CUDAallocator_device_mem_allocated;

inline size_t getCUDAdeviceMemAllocated() { return CUDAallocator_device_mem_allocated; }

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

template<class T1, class T2>
bool operator==(const CUDAManagedAllocator<T1>&, const CUDAManagedAllocator<T2>&)
{
  return true;
}
template<class T1, class T2>
bool operator!=(const CUDAManagedAllocator<T1>&, const CUDAManagedAllocator<T2>&)
{
  return false;
}


/** allocator for CUDA device memory
 * @tparam T data type
 *
 * using this with something other than Ohmms containers?
 *  -- use caution, write unit tests! --
 * It's not tested beyond use in some unit tests using std::vector with constant size.
 * CUDAAllocator appears to meet all the nonoptional requirements of a c++ Allocator.
 *
 * Some of the default implementations in std::allocator_traits
 * of optional Allocator requirements may cause runtime or compilation failures.
 * They assume there is only one memory space and that the host has access to it.
 */
template<typename T>
class CUDAAllocator
{
public:
  using value_type    = T;
  using size_type     = size_t;
  using pointer       = T*;
  using const_pointer = const T*;

  CUDAAllocator() = default;
  template<class U>
  CUDAAllocator(const CUDAAllocator<U>&)
  {}

  template<class U>
  struct rebind
  {
    using other = CUDAAllocator<U>;
  };

  T* allocate(std::size_t n)
  {
    void* pt;
    cudaErrorCheck(cudaMalloc(&pt, n * sizeof(T)), "Allocation failed in CUDAAllocator!");
    CUDAallocator_device_mem_allocated += n * sizeof(T);
    return static_cast<T*>(pt);
  }
  void deallocate(T* p, std::size_t n)
  {
    cudaErrorCheck(cudaFree(p), "Deallocation failed in CUDAAllocator!");
    CUDAallocator_device_mem_allocated -= n * sizeof(T);
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
    cudaErrorCheck(cudaMemcpy(device_ptr, host_ptr, sizeof(T) * n, cudaMemcpyHostToDevice),
                   "cudaMemcpy failed in copyToDevice");
  }

  void copyFromDevice(T* host_ptr, T* device_ptr, size_t n)
  {
    cudaErrorCheck(cudaMemcpy(host_ptr, device_ptr, sizeof(T) * n, cudaMemcpyDeviceToHost),
                   "cudaMemcpy failed in copyFromDevice");
  }

  void copyDeviceToDevice(T* to_ptr, size_t n, T* from_ptr)
  {
    cudaErrorCheck(cudaMemcpy(to_ptr, from_ptr, sizeof(T) * n, cudaMemcpyDeviceToDevice),
                   "cudaMemcpy failed in copyDeviceToDevice");
  }
};

template<class T1, class T2>
bool operator==(const CUDAAllocator<T1>&, const CUDAAllocator<T2>&)
{
  return true;
}
template<class T1, class T2>
bool operator!=(const CUDAAllocator<T1>&, const CUDAAllocator<T2>&)
{
  return false;
}
template<typename T>

struct qmc_allocator_traits<qmcplusplus::CUDAAllocator<T>>
{
  static const bool is_host_accessible = false;
  static const bool is_dual_space      = false;
  static void fill_n(T* ptr, size_t n, const T& value) { qmcplusplus::CUDAfill_n(ptr, n, value); }
};

/** allocator for CUDA host pinned memory
 * @tparam T data type
 */
template<typename T>
struct CUDAHostAllocator
{
  using value_type    = T;
  using size_type     = size_t;
  using pointer       = T*;
  using const_pointer = const T*;

  CUDAHostAllocator() = default;
  template<class U>
  CUDAHostAllocator(const CUDAHostAllocator<U>&)
  {}

  template<class U>
  struct rebind
  {
    using other = CUDAHostAllocator<U>;
  };

  T* allocate(std::size_t n)
  {
    void* pt;
    cudaErrorCheck(cudaMallocHost(&pt, n * sizeof(T)), "Allocation failed in CUDAHostAllocator!");
    return static_cast<T*>(pt);
  }
  void deallocate(T* p, std::size_t) { cudaErrorCheck(cudaFreeHost(p), "Deallocation failed in CUDAHostAllocator!"); }
};

template<class T1, class T2>
bool operator==(const CUDAHostAllocator<T1>&, const CUDAHostAllocator<T2>&)
{
  return true;
}
template<class T1, class T2>
bool operator!=(const CUDAHostAllocator<T1>&, const CUDAHostAllocator<T2>&)
{
  return false;
}

/** allocator locks memory pages allocated by ULPHA
 * @tparam T data type
 * @tparam ULPHA host memory allocator using unlocked page
 *
 * ULPHA cannot be CUDAHostAllocator
 */
template<typename T, class ULPHA = std::allocator<T>>
struct CUDALockedPageAllocator : public ULPHA
{
  using value_type    = typename ULPHA::value_type;
  using size_type     = typename ULPHA::size_type;
  using pointer       = typename ULPHA::pointer;
  using const_pointer = typename ULPHA::const_pointer;

  CUDALockedPageAllocator() = default;
  template<class U, class V>
  CUDALockedPageAllocator(const CUDALockedPageAllocator<U, V>&)
  {}

  template<class U, class V>
  struct rebind
  {
    using other = CUDALockedPageAllocator<U, V>;
  };

  value_type* allocate(std::size_t n)
  {
    static_assert(std::is_same<T, value_type>::value, "CUDALockedPageAllocator and ULPHA data types must agree!");
    value_type* pt = ULPHA::allocate(n);
    cudaErrorCheck(cudaHostRegister(pt, n * sizeof(T), cudaHostRegisterDefault),
                   "cudaHostRegister failed in CUDALockedPageAllocator!");
    return pt;
  }

  void deallocate(value_type* pt, std::size_t n)
  {
    cudaErrorCheck(cudaHostUnregister(pt), "cudaHostUnregister failed in CUDALockedPageAllocator!");
    ULPHA::deallocate(pt, n);
  }
};

} // namespace qmcplusplus

#endif
