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
/** @file SYCLallocator.hpp
 * this file provides three C++ memory allocators using SYCL specific memory allocation functions.
 *
 * SYCLManagedAllocator allocates SYCL unified memory
 * SYCLAllocator allocates SYCL device memory
 * SYCLHostAllocator allocates SYCL host pinned memory
 */
#ifndef QMCPLUSPLUS_SYCL_ALLOCATOR_H
#define QMCPLUSPLUS_SYCL_ALLOCATOR_H

#include <stdexcept>
#include "allocator_traits.hpp"
#include "CPU/SIMD/alignment.config.h"
#include <CL/sycl.hpp>
namespace sycl = cl::sycl;

namespace qmcplusplus
{
struct SYCLAllocatorBase
{
public:
  // ensure same context is used for all allocators/types and kernel execution
  // the quick and dirty implentation for now is to just use one queue for everything
  // eventually, this should be changed to a single shared context allowing multiple queues
  inline static sycl::queue q;
};
/** allocator for SYCL unified memory
 * @tparm T data type
 */
template<typename T>
struct SYCLManagedAllocator : SYCLAllocatorBase
{
  typedef T value_type;
  typedef size_t size_type;
  typedef T* pointer;
  typedef const T* const_pointer;

  SYCLManagedAllocator() = default;
  template<class U>
  SYCLManagedAllocator(const SYCLManagedAllocator<U>&)
  {}

  template<class U>
  struct rebind
  {
    typedef SYCLManagedAllocator<U> other;
  };

  T* allocate(std::size_t n)
  {
    sycl::device dev  = q.get_device();
    sycl::context ctx = q.get_context();
    return static_cast<T*>(sycl::malloc_shared(n * sizeof(T), dev, ctx));
  }
  void deallocate(T* p, std::size_t)
  {
    sycl::device dev  = q.get_device();
    sycl::context ctx = q.get_context();
    sycl::free(p, ctx);
  }
};

template<class T1, class T2>
bool operator==(const SYCLManagedAllocator<T1>&, const SYCLManagedAllocator<T2>&)
{
  return true;
}
template<class T1, class T2>
bool operator!=(const SYCLManagedAllocator<T1>&, const SYCLManagedAllocator<T2>&)
{
  return false;
}

/** allocator for SYCL device memory
 * @tparm T data type
 */
template<typename T>
struct SYCLAllocator
{
  typedef T value_type;
  typedef size_t size_type;
  typedef T* pointer;
  typedef const T* const_pointer;

  SYCLAllocator() = default;
  template<class U>
  SYCLAllocator(const SYCLAllocator<U>&)
  {}

  template<class U>
  struct rebind
  {
    typedef SYCLAllocator<U> other;
  };

  T* allocate(std::size_t n) { throw("SYCLAllocator::allocate not yet implemented"); }
  void deallocate(T* p, std::size_t) { throw("SYCLAllocator::deallocate not yet implemented"); }
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
struct allocator_traits<SYCLAllocator<T>>
{
  const static bool is_host_accessible = false;
};

/** allocator for SYCL host pinned memory
 * @tparm T data type
 */
template<typename T>
struct SYCLHostAllocator
{
  typedef T value_type;
  typedef size_t size_type;
  typedef T* pointer;
  typedef const T* const_pointer;

  SYCLHostAllocator() = default;
  template<class U>
  SYCLHostAllocator(const SYCLHostAllocator<U>&)
  {}

  template<class U>
  struct rebind
  {
    typedef SYCLHostAllocator<U> other;
  };

  T* allocate(std::size_t n) { throw("SYCLHostAllocator::allocate not yet implemented"); }
  void deallocate(T* p, std::size_t) { throw("SYCLHostAllocator::deallocate not yet implemented"); }
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

/** allocator locks memory pages allocated by ULPHA
 * @tparm T data type
 * @tparm ULPHA host memory allocator using unlocked page
 *
 * ULPHA cannot be SYCLHostAllocator
 */
template<typename T, class ULPHA = std::allocator<T>>
struct SYCLLockedPageAllocator : public ULPHA
{
  using value_type    = typename ULPHA::value_type;
  using size_type     = typename ULPHA::size_type;
  using pointer       = typename ULPHA::pointer;
  using const_pointer = typename ULPHA::const_pointer;

  SYCLLockedPageAllocator() = default;
  template<class U, class V>
  SYCLLockedPageAllocator(const SYCLLockedPageAllocator<U, V>&)
  {}
  template<class U, class V>
  struct rebind
  {
    typedef SYCLLockedPageAllocator<U, V> other;
  };

  value_type* allocate(std::size_t n) { throw("SYCLLockedPageAllocator::allocate not yet implemented"); }

  void deallocate(value_type* pt, std::size_t n) { throw("SYCLLockedPageAllocator::dellalocate not yet implemented"); }
};

} // namespace qmcplusplus

#endif
