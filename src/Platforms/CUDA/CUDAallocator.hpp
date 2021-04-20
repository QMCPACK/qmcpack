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
#include <cuda_runtime_api.h>
#include "config.h"
#include "cudaError.h"
#include "allocator_traits.hpp"
#include "CUDAfill.hpp"

namespace qmcplusplus
{
extern std::atomic<size_t> CUDAallocator_device_mem_allocated;

inline size_t getCUDAdeviceMemAllocated() { return CUDAallocator_device_mem_allocated; }

/** allocator for CUDA unified memory
 * @tparm T data type
 */
template<typename T>
struct CUDAManagedAllocator
{
  typedef T value_type;
  typedef size_t size_type;
  typedef T* pointer;
  typedef const T* const_pointer;

  CUDAManagedAllocator() = default;
  template<class U>
  CUDAManagedAllocator(const CUDAManagedAllocator<U>&)
  {}

  template<class U>
  struct rebind
  {
    typedef CUDAManagedAllocator<U> other;
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
 * @tparm T data type
 *
 *  using this with something other than Ohmms containers?
 *   -- use caution --
 *  it does not provide all the allocator_traits specialization and methods for an
 *  std::allocator that some std containers expect.
 *  This can result in unexpected behavior at runtime or if you are lucky compilation errors.
 */
template<typename T>
struct CUDAAllocator
{
  typedef T value_type;
  typedef size_t size_type;
  typedef T* pointer;
  typedef const T* const_pointer;

  CUDAAllocator() = default;
  template<class U>
  CUDAAllocator(const CUDAAllocator<U>&)
  {}

  template<class U>
  struct rebind
  {
    typedef CUDAAllocator<U> other;
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

  /** Don't do anything on construct, we don't initialize or otherwise touch device memory at construction.
   *  std::vector at least wants to call the default initializer on each element using this.
   *
   *  The standard is a bit confusing on this point. This is an optional requirement of Allocator from C++11 on.
   *  Its deprecated for the std::allocator in c++17 and will be removed in c++20.  We are not implementing
   *  std::allocator.
   *
   *  Since std::allocator_traits handles the case where no construct method is present in the Allocator.
   *  But for allocators other than std::allocator the std::allocator_traits will call the Allocators construct
   *  method.
   */
  template<class U, class... Args>
  static void construct(U* p, Args&&... args)
  {}

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
 * @tparm T data type
 */
template<typename T>
struct CUDAHostAllocator
{
  typedef T value_type;
  typedef size_t size_type;
  typedef T* pointer;
  typedef const T* const_pointer;

  CUDAHostAllocator() = default;
  template<class U>
  CUDAHostAllocator(const CUDAHostAllocator<U>&)
  {}

  template<class U>
  struct rebind
  {
    typedef CUDAHostAllocator<U> other;
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
 * @tparm T data type
 * @tparm ULPHA host memory allocator using unlocked page
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
    typedef CUDALockedPageAllocator<U, V> other;
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
