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

#ifndef QMCPLUSPLUS_MEMMANAGE_H
#define QMCPLUSPLUS_MEMMANAGE_H

#include "PlatformKinds.hpp"

namespace qmcplusplus
{

namespace compute
{

template<PlatformKind PL>
class MemManage;

/** allocator for device memory
 * @tparam T data type
 *
 * using this with something other than Ohmms containers?
 *  -- use caution, write unit tests! --
 * It's not tested beyond use in some unit tests using std::vector with constant size.
 * DeviceAllocatorImpl appears to meet all the nonoptional requirements of a c++ Allocator.
 *
 * Some of the default implementations in std::allocator_traits
 * of optional Allocator requirements may cause runtime or compilation failures.
 * They assume there is only one memory space and that the host has access to it.
 */
template<typename T, PlatformKind PL>
class DeviceAllocatorImpl
{
public:
  using value_type    = T;
  using size_type     = size_t;
  using pointer       = T*;
  using const_pointer = const T*;

  DeviceAllocatorImpl() = default;
  template<class U>
  DeviceAllocatorImpl(const DeviceAllocatorImpl<U, PL>&)
  {}

  template<class U>
  struct rebind
  {
    using other = DeviceAllocatorImpl<U, PL>;
  };

  T* allocate(std::size_t n)
  {
    void* pt;
    MemManage<PL>::mallocDevice(&pt, n * sizeof(T));
    MemManage<PL>::device_mem_allocated_ += n * sizeof(T);
    return static_cast<T*>(pt);
  }

  void deallocate(T* p, std::size_t n)
  {
    MemManage<PL>::freeDevice(p);
    MemManage<PL>::device_mem_allocated_ -= n * sizeof(T);
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

  static void memcpy(T* dst, const T* src, size_t size) { MemManage<PL>::memcpy(dst, src, sizeof(T) * size); }
};

/** allocator for host memory
 * @tparam T data type
 */
template<typename T, PlatformKind PL>
struct HostAllocatorImpl
{
  using value_type    = T;
  using size_type     = size_t;
  using pointer       = T*;
  using const_pointer = const T*;

  HostAllocatorImpl() = default;
  template<class U>
  HostAllocatorImpl(const HostAllocatorImpl<U, PL>&)
  {}

  template<class U>
  struct rebind
  {
    using other = HostAllocatorImpl<U, PL>;
  };

  T* allocate(std::size_t n)
  {
    void* pt;
    MemManage<PL>::mallocHost(&pt, n * sizeof(T));
    return static_cast<T*>(pt);
  }

  void deallocate(T* p, std::size_t) { MemManage<PL>::freeHost(p); }
};

/** allocator locks memory pages allocated by ULPHA
 * @tparam T data type
 * @tparam ULPHA host memory allocator using unlocked page
 *
 * ULPHA cannot be HostAllocator
 */
template<typename T, PlatformKind PL, class ULPHA = std::allocator<T>>
struct PageLockedAllocatorImpl : public ULPHA
{
  using value_type    = typename ULPHA::value_type;
  using size_type     = typename ULPHA::size_type;
  using pointer       = typename ULPHA::pointer;
  using const_pointer = typename ULPHA::const_pointer;

  PageLockedAllocatorImpl() = default;
  template<class U, class V>
  PageLockedAllocatorImpl(const PageLockedAllocatorImpl<U, PL, V>&)
  {}

  template<class U>
  struct rebind
  {
    using other = PageLockedAllocatorImpl<U, PL, typename std::allocator_traits<ULPHA>::template rebind_alloc<U>>;
  };

  value_type* allocate(std::size_t n)
  {
    static_assert(std::is_same<T, value_type>::value, "PageLockedAllocatorImpl and ULPHA data types must agree!");
    value_type* pt = ULPHA::allocate(n);
    MemManage<PL>::registerHost(pt, n * sizeof(T));
    return pt;
  }

  void deallocate(value_type* pt, std::size_t n)
  {
    MemManage<PL>::unregisterHost(pt);
    ULPHA::deallocate(pt, n);
  }
};
} // namespace compute

} // namespace qmcplusplus

#endif
