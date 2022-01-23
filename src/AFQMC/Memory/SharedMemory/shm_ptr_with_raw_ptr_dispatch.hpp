//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//    Lawrence Livermore National Laboratory
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_SHM_PTR_WITH_RAW_PTR_DISPATCH_HPP
#define AFQMC_SHM_PTR_WITH_RAW_PTR_DISPATCH_HPP

#include <functional>
#include "Configuration.h"
#include "AFQMC/config.0.h"
#include <cassert>

#include "multi/array_ref.hpp"
#include "mpi3/communicator.hpp"
#include "mpi3/group.hpp"
#include "mpi3/shared_communicator.hpp"
#include "mpi3/shared_window.hpp"
#include "AFQMC/Memory/raw_pointers.hpp"

namespace shm
{
using qmcplusplus::afqmc::to_address;

namespace mpi3  = boost::mpi3;
namespace multi = boost::multi;

template<class T>
struct allocator_shm_ptr_with_raw_ptr_dispatch;
template<class T>
struct shm_ptr_with_raw_ptr_dispatch;
template<typename T, typename Size>
shm_ptr_with_raw_ptr_dispatch<T> destroy_n(shm_ptr_with_raw_ptr_dispatch<T> first, Size n);
template<typename T, typename Size, typename TT>
shm_ptr_with_raw_ptr_dispatch<T> uninitialized_fill_n(shm_ptr_with_raw_ptr_dispatch<T> first, Size n, TT const& val);

template<>
struct shm_ptr_with_raw_ptr_dispatch<const void>
{
  using T = const void;
  std::shared_ptr<mpi3::shared_window<char>> wSP_;
  std::ptrdiff_t offset = 0;
  shm_ptr_with_raw_ptr_dispatch(std::nullptr_t = nullptr) {}
  shm_ptr_with_raw_ptr_dispatch(shm_ptr_with_raw_ptr_dispatch const& other)            = default;
  shm_ptr_with_raw_ptr_dispatch& operator=(shm_ptr_with_raw_ptr_dispatch const& other) = default;
  shm_ptr_with_raw_ptr_dispatch& operator=(std::nullptr_t)
  {
    wSP_.reset();
    return *this;
  }
  bool operator==(std::nullptr_t) const { return (bool)wSP_; }
  bool operator!=(std::nullptr_t) const { return not operator==(nullptr); }

private:
  shm_ptr_with_raw_ptr_dispatch(std::shared_ptr<mpi3::shared_window<char>> wSP) : wSP_{wSP} {}
  template<class>
  friend struct shm_ptr_with_raw_ptr_dispatch;
};

template<>
struct shm_ptr_with_raw_ptr_dispatch<void>
{
  using T            = void;
  using element_type = T;
  std::shared_ptr<mpi3::shared_window<char>> wSP_;
  std::ptrdiff_t offset = 0;
  shm_ptr_with_raw_ptr_dispatch(std::nullptr_t = nullptr) {}
  shm_ptr_with_raw_ptr_dispatch(shm_ptr_with_raw_ptr_dispatch const& other)            = default;
  shm_ptr_with_raw_ptr_dispatch& operator=(shm_ptr_with_raw_ptr_dispatch const& other) = default;
  shm_ptr_with_raw_ptr_dispatch& operator=(std::nullptr_t)
  {
    wSP_.reset();
    return *this;
  }
  template<typename Q>
  shm_ptr_with_raw_ptr_dispatch(shm_ptr_with_raw_ptr_dispatch<Q> const& ptr) : wSP_(ptr.wSP_), offset(ptr.offset)
  {}
  bool operator==(std::nullptr_t) const { return (bool)wSP_; }
  bool operator!=(std::nullptr_t) const { return not operator==(nullptr); }
  explicit operator bool() const { return (bool)wSP_; }
  bool operator==(shm_ptr_with_raw_ptr_dispatch const& other) const
  {
    return (wSP_->base(0) == other.wSP_->base(0)) && (offset == other.offset);
  }
  bool operator!=(shm_ptr_with_raw_ptr_dispatch const& other) const { return not(*this == other); }
  bool operator<(shm_ptr_with_raw_ptr_dispatch const& other) const
  {
    return wSP_->base(0) + offset < other.wSP_->base(0) + other.offset;
  }
};

template<class T>
struct allocator_shm_ptr_with_raw_ptr_dispatch;

template<class T>
struct shm_ptr_with_raw_ptr_dispatch
{
  using element_type           = T;
  using difference_type        = std::ptrdiff_t;
  using value_type             = std::decay_t<T>; // std::remove_cv_t<T>; // T until C++20?
  using pointer                = T*;              // TODO self?
  using reference              = T&;              //TODO fancy_reference?
  using iterator_category      = std::random_access_iterator_tag;
  using rebind_const           = shm_ptr_with_raw_ptr_dispatch<const T>;
  using default_allocator_type = allocator_shm_ptr_with_raw_ptr_dispatch<value_type>;
  std::shared_ptr<mpi3::shared_window<char>> wSP_;
  std::ptrdiff_t offset = 0; // in Bytes
  shm_ptr_with_raw_ptr_dispatch() {}
  shm_ptr_with_raw_ptr_dispatch(std::nullptr_t) {}
  shm_ptr_with_raw_ptr_dispatch(shm_ptr_with_raw_ptr_dispatch const& other) = default;
  template<typename Q>
  shm_ptr_with_raw_ptr_dispatch(shm_ptr_with_raw_ptr_dispatch<Q> const& ptr) : wSP_(ptr.wSP_), offset(ptr.offset)
  {}
  shm_ptr_with_raw_ptr_dispatch& operator=(shm_ptr_with_raw_ptr_dispatch const& other) = default;
  shm_ptr_with_raw_ptr_dispatch& operator=(std::nullptr_t) { return *this; }
  ~shm_ptr_with_raw_ptr_dispatch() = default;
  T& operator*() const { return *(reinterpret_cast<T*>(wSP_->base(0) + offset)); }
  T& operator[](int idx) const { return (reinterpret_cast<T*>(wSP_->base(0) + offset))[idx]; }
  T* operator->() const { return reinterpret_cast<T*>(wSP_->base(0) + offset); }
  operator shm_ptr_with_raw_ptr_dispatch<const void>() const { return shm_ptr_with_raw_ptr_dispatch<const void>{wSP_}; }
  T* get() const
  {
    if (wSP_ == nullptr)
      return nullptr;
    else
      return reinterpret_cast<T*>(wSP_->base(0) + offset);
  }
  explicit operator T*() const { return get(); }
  explicit operator bool() const { return (bool)wSP_; } //.get();}
  bool operator==(std::nullptr_t) const { return not(bool) wSP_; }
  bool operator!=(std::nullptr_t) const { return not operator==(nullptr); }
  /*
        operator shm_ptr_with_raw_ptr_dispatch<T const>() const{
                shm_ptr_with_raw_ptr_dispatch<T const> ret;
                ret.wSP_ = wSP_;
                ret.offset = offset;
                return ret;
        }
        operator shm_ptr_with_raw_ptr_dispatch<void const>() const{
                shm_ptr_with_raw_ptr_dispatch<void const> ret;
                ret.wSP_ = wSP_;
                ret.offset = offset;
                return ret;
        }
*/
  shm_ptr_with_raw_ptr_dispatch operator+(std::ptrdiff_t d) const
  {
    shm_ptr_with_raw_ptr_dispatch ret(*this);
    ret += d;
    return ret;
  }
  shm_ptr_with_raw_ptr_dispatch operator-(std::ptrdiff_t d) const
  {
    shm_ptr_with_raw_ptr_dispatch ret(*this);
    ret -= d;
    return ret;
  }
  std::ptrdiff_t operator-(shm_ptr_with_raw_ptr_dispatch<T> other) const
  {
    std::ptrdiff_t pdiff(wSP_->base(0) + offset - (other.wSP_->base(0) + other.offset));
    if (pdiff % sizeof(T) != 0)
    {
      throw;
    }
    return pdiff / sizeof(T);
  }
  shm_ptr_with_raw_ptr_dispatch& operator--()
  {
    offset -= sizeof(T);
    return *this;
  }
  shm_ptr_with_raw_ptr_dispatch& operator++()
  {
    offset += sizeof(T);
    return *this;
  }
  shm_ptr_with_raw_ptr_dispatch& operator-=(std::ptrdiff_t d)
  {
    offset -= d * sizeof(T);
    return *this;
  }
  shm_ptr_with_raw_ptr_dispatch& operator+=(std::ptrdiff_t d)
  {
    offset += d * sizeof(T);
    return *this;
  }
  bool operator==(shm_ptr_with_raw_ptr_dispatch<T> const& other) const
  {
    return wSP_->base(0) == other.wSP_->base(0) and offset == other.offset;
  }
  bool operator!=(shm_ptr_with_raw_ptr_dispatch<T> const& other) const { return not((*this) == other); }
  bool operator<(shm_ptr_with_raw_ptr_dispatch<T> const& other) const
  {
    return wSP_->base(0) + offset < other.wSP_->base(0) + other.offset;
  }
  static element_type* to_address(shm_ptr_with_raw_ptr_dispatch p) noexcept
  {
    return reinterpret_cast<element_type*>(p.wSP_->base(0) + p.offset);
  }
  friend pointer to_address(shm_ptr_with_raw_ptr_dispatch const& p)
  {
    return shm_ptr_with_raw_ptr_dispatch::to_address(p);
  }
};

template<class T = void>
struct allocator_shm_ptr_with_raw_ptr_dispatch
{
  template<class U>
  struct rebind
  {
    using other = allocator_shm_ptr_with_raw_ptr_dispatch<U>;
  };
  using value_type      = T;
  using pointer         = shm_ptr_with_raw_ptr_dispatch<T>;
  using const_pointer   = shm_ptr_with_raw_ptr_dispatch<T const>;
  using size_type       = mpi3::size_t;                  // std::size_t;
  using difference_type = std::make_signed_t<size_type>; //std::ptrdiff_t;

  mpi3::shared_communicator* commP_;

  allocator_shm_ptr_with_raw_ptr_dispatch() = delete;
  allocator_shm_ptr_with_raw_ptr_dispatch(mpi3::shared_communicator& comm) : commP_(&comm) {}
  allocator_shm_ptr_with_raw_ptr_dispatch(allocator_shm_ptr_with_raw_ptr_dispatch const& other) : commP_(other.commP_)
  {}
  ~allocator_shm_ptr_with_raw_ptr_dispatch() = default;
  template<class U>
  allocator_shm_ptr_with_raw_ptr_dispatch(allocator_shm_ptr_with_raw_ptr_dispatch<U> const& o) : commP_(o.commP_)
  {}

  shm_ptr_with_raw_ptr_dispatch<T> allocate(size_type n, const void* /*hint*/ = 0)
  {
    shm_ptr_with_raw_ptr_dispatch<T> ret = 0;
    ret.wSP_.reset(new mpi3::shared_window<char>{*commP_, commP_->root() ? (long(n * sizeof(T))) : 0, int(sizeof(T))});
    return ret;
  }
  void deallocate(shm_ptr_with_raw_ptr_dispatch<T> ptr, size_type) { ptr.wSP_.reset(); }
  allocator_shm_ptr_with_raw_ptr_dispatch& operator=(allocator_shm_ptr_with_raw_ptr_dispatch const& other)
  {
    assert((*this) == other); // TODO make comm a shared_ptr
    return *this;
  }
  bool operator==(allocator_shm_ptr_with_raw_ptr_dispatch const& other) const { return commP_ == other.commP_; }
  bool operator!=(allocator_shm_ptr_with_raw_ptr_dispatch const& other) const { return not(other == *this); }
  // this routine synchronizes
  template<class U, class... As>
  void construct(U p, As&&... as)
  {
    uninitialized_fill_n(p, 1, std::forward<As>(as)...);
  }
  // this routine synchronizes
  template<class U, class... As>
  void destroy(U p, As&&... as)
  {
    destroy_n(p, 1);
  }
};

struct memory_resource_shm_ptr_with_raw_ptr_dispatch
{
  using pointer = shm_ptr_with_raw_ptr_dispatch<void>;

  mpi3::shared_communicator* commP_;

  memory_resource_shm_ptr_with_raw_ptr_dispatch(memory_resource_shm_ptr_with_raw_ptr_dispatch const&) = default;

  shm_ptr_with_raw_ptr_dispatch<void> allocate(std::size_t size, std::size_t alignment = alignof(std::max_align_t))
  {
    shm_ptr_with_raw_ptr_dispatch<char> ret = 0;
    ret.wSP_.reset(new mpi3::shared_window<char>{*commP_, commP_->root() ? long(size) : 0, int(alignment)});
    return ret;
  }
  void deallocate(shm_ptr_with_raw_ptr_dispatch<void> ptr, std::size_t) { ptr.wSP_.reset(); }

  bool operator==(memory_resource_shm_ptr_with_raw_ptr_dispatch const& other) const { return commP_ == other.commP_; }

  bool operator!=(memory_resource_shm_ptr_with_raw_ptr_dispatch const& other) const { return not(other == *this); }
};

template<typename T, typename Size, typename TT>
shm_ptr_with_raw_ptr_dispatch<T> fill_n(shm_ptr_with_raw_ptr_dispatch<T> first, Size n, TT const& val)
{
  if (n == 0)
    return first;
  first.wSP_->fence();
  if (first.wSP_->get_group().root())
    std::fill_n(to_address(first), n, val);
  first.wSP_->fence();
  first.wSP_->fence();
  mpi3::communicator(first.wSP_->get_group(), 0).barrier();
  return first + n;
}

template<typename T, typename Q>
shm_ptr_with_raw_ptr_dispatch<T> fill(shm_ptr_with_raw_ptr_dispatch<T> first,
                                      shm_ptr_with_raw_ptr_dispatch<T> last,
                                      Q const& val)
{
  if (std::distance(first, last) == 0)
    return first;
  return fill_n(first, std::distance(first, last), val);
}

template<typename T, typename Size, typename TT>
shm_ptr_with_raw_ptr_dispatch<T> uninitialized_fill_n(shm_ptr_with_raw_ptr_dispatch<T> first, Size n, TT const& val)
{
  if (n == 0)
    return first;
  first.wSP_->fence();
  if (first.wSP_->get_group().root())
    std::uninitialized_fill_n(to_address(first), n, val);
  first.wSP_->fence();
  first.wSP_->fence();
  mpi3::communicator(first.wSP_->get_group(), 0).barrier();
  return first + n;
}

template<class Alloc, typename T, typename Size, typename TT>
shm_ptr_with_raw_ptr_dispatch<T> uninitialized_fill_n(Alloc& a,
                                                      shm_ptr_with_raw_ptr_dispatch<T> first,
                                                      Size n,
                                                      TT const& val)
{
  if (n == 0)
    return first;
  first.wSP_->fence();
  if (first.wSP_->get_group().root())
    std::uninitialized_fill_n(to_address(first), n, val); // change to to_pointer
  first.wSP_->fence();
  first.wSP_->fence();
  mpi3::communicator(first.wSP_->get_group(), 0).barrier();
  return first + n;
}

template<class Alloc, typename T, typename Size, typename TT>
shm_ptr_with_raw_ptr_dispatch<T> alloc_uninitialized_fill_n(Alloc& a,
                                                            shm_ptr_with_raw_ptr_dispatch<T> first,
                                                            Size n,
                                                            TT const& val)
{
  if (n == 0)
    return first;
  first.wSP_->fence();
  if (first.wSP_->get_group().root())
    std::uninitialized_fill_n(to_address(first), n, val); // change to to_pointer
  first.wSP_->fence();
  first.wSP_->fence();
  mpi3::communicator(first.wSP_->get_group(), 0).barrier();
  return first + n;
}

template<class Alloc, typename T, typename Size, typename TT>
shm_ptr_with_raw_ptr_dispatch<T> alloc_uninitialized_fill(Alloc& a,
                                                          shm_ptr_with_raw_ptr_dispatch<T> first,
                                                          shm_ptr_with_raw_ptr_dispatch<T> last,
                                                          TT const& val)
{
  if (std::distance(first, last) == 0)
    return first;
  return alloc_uninitialized_fill_n(a, first, std::distance(first, last), val);
}

template<typename T, typename Size>
shm_ptr_with_raw_ptr_dispatch<T> destroy_n(shm_ptr_with_raw_ptr_dispatch<T> first, Size n)
{
  if (n == 0)
    return first;
  first.wSP_->fence();
  if (first.wSP_->get_group().root())
  {
    auto first_ptr = to_address(first);
    for (; n > 0; (void)++first_ptr, --n)
      first_ptr->~T();
  }
  first.wSP_->fence();
  first.wSP_->fence();
  mpi3::communicator(first.wSP_->get_group(), 0).barrier();
  return first + n;
}

template<class Alloc, typename T, typename Size>
shm_ptr_with_raw_ptr_dispatch<T> destroy_n(Alloc& a, shm_ptr_with_raw_ptr_dispatch<T> first, Size n)
{
  if (n == 0)
    return first;
  first.wSP_->fence();
  if (first.wSP_->get_group().root())
  {
    auto first_ptr = to_address(first);
    for (; n > 0; (void)++first_ptr, --n)
      first_ptr->~T();
  }
  first.wSP_->fence();
  first.wSP_->fence();
  mpi3::communicator(first.wSP_->get_group(), 0).barrier();
  return first + n;
}

template<class Alloc, typename T, typename Size>
shm_ptr_with_raw_ptr_dispatch<T> alloc_destroy_n(Alloc& a, shm_ptr_with_raw_ptr_dispatch<T> first, Size n)
{
  if (n == 0)
    return first;
  first.wSP_->fence();
  if (first.wSP_->get_group().root())
  {
    auto first_ptr = to_address(first);
    for (; n > 0; (void)++first_ptr, --n)
      first_ptr->~T();
  }
  first.wSP_->fence();
  first.wSP_->fence();
  mpi3::communicator(first.wSP_->get_group(), 0).barrier();
  return first + n;
}

template<class It1, typename T, typename Size>
shm_ptr_with_raw_ptr_dispatch<T> copy_n_cast(It1 first, Size n, shm_ptr_with_raw_ptr_dispatch<T> d_first)
{
  if (n == 0)
    return d_first;
  d_first.wSP_->fence();
  using qmcplusplus::afqmc::copy_n_cast;
  if (d_first.wSP_->get_group().root())
    copy_n_cast(to_address(first), n, to_address(d_first));
  d_first.wSP_->fence();
  mpi3::communicator(d_first.wSP_->get_group(), 0).barrier();
  return d_first + n;
}

template<class It1, typename T, typename Size>
shm_ptr_with_raw_ptr_dispatch<T> copy_n(It1 first, Size n, shm_ptr_with_raw_ptr_dispatch<T> d_first)
{
  if (n == 0)
    return d_first;
  d_first.wSP_->fence();
  using std::copy_n;
  if (d_first.wSP_->get_group().root())
    copy_n(first, n, to_address(d_first));
  d_first.wSP_->fence();
  mpi3::communicator(d_first.wSP_->get_group(), 0).barrier();
  return d_first + n;
}

template<class It1, typename T>
shm_ptr_with_raw_ptr_dispatch<T> copy(It1 first, It1 last, shm_ptr_with_raw_ptr_dispatch<T> d_first)
{
  using std::distance;
  if (distance(first, last) == 0)
    return d_first;
  first.wSP_->fence();
  using std::copy;
  if (d_first.wSP_->get_group().root())
    copy(first, last, to_address(d_first));
  first.wSP_->fence();
  mpi3::communicator(d_first.wSP_->get_group(), 0).barrier();
  return d_first + distance(first, last);
}

template<class It1, class Size, typename T>
shm_ptr_with_raw_ptr_dispatch<T> uninitialized_copy_n(It1 f, Size n, shm_ptr_with_raw_ptr_dispatch<T> d)
{
  if (n == 0)
    return d;
  d.wSP_->fence();
  using std::uninitialized_copy_n;
  if (d.wSP_->get_group().root())
    uninitialized_copy_n(f, n, to_address(d));
  d.wSP_->fence();
  mpi3::communicator(d.wSP_->get_group(), 0).barrier();
  return d + n;
}

template<class Alloc, class It1, class Size, typename T>
shm_ptr_with_raw_ptr_dispatch<T> uninitialized_copy_n(Alloc& a, It1 f, Size n, shm_ptr_with_raw_ptr_dispatch<T> d)
{
  if (n == 0)
    return d;
  d.wSP_->fence();
  using std::uninitialized_copy_n;
  if (d.wSP_->get_group().root())
    uninitialized_copy_n(f, n, to_address(d));
  d.wSP_->fence();
  mpi3::communicator(d.wSP_->get_group(), 0).barrier();
  return d + n;
}

template<class Alloc, class It1, class Size, typename T>
shm_ptr_with_raw_ptr_dispatch<T> alloc_uninitialized_copy_n(Alloc& a, It1 f, Size n, shm_ptr_with_raw_ptr_dispatch<T> d)
{
  if (n == 0)
    return d;
  d.wSP_->fence();
  using std::uninitialized_copy_n;
  if (d.wSP_->get_group().root())
    uninitialized_copy_n(f, n, to_address(d));
  d.wSP_->fence();
  mpi3::communicator(d.wSP_->get_group(), 0).barrier();
  return d + n;
}

template<class It1, typename T>
shm_ptr_with_raw_ptr_dispatch<T> uninitialized_copy(It1 f, It1 l, shm_ptr_with_raw_ptr_dispatch<T> d)
{
  using std::distance;
  if (distance(f, l) == 0)
    return d;
  d.wSP_->fence();
  using std::uninitialized_copy;
  if (d.wSP_->get_group().root())
    uninitialized_copy(f, l, to_address(d));
  d.wSP_->fence();
  mpi3::communicator(d.wSP_->get_group(), 0).barrier();
  return d + distance(f, l);
}

template<class Alloc, class It1, typename T>
shm_ptr_with_raw_ptr_dispatch<T> uninitialized_copy(Alloc& a, It1 f, It1 l, shm_ptr_with_raw_ptr_dispatch<T> d)
{
  using std::distance;
  if (distance(f, l) == 0)
    return d;
  d.wSP_->fence();
  using std::uninitialized_copy;
  if (d.wSP_->get_group().root())
    uninitialized_copy(f, l, to_address(d));
  d.wSP_->fence();
  mpi3::communicator(d.wSP_->get_group(), 0).barrier();
  return d + distance(f, l);
}

template<class Alloc, class It1, typename T>
shm_ptr_with_raw_ptr_dispatch<T> alloc_uninitialized_copy(Alloc& a, It1 f, It1 l, shm_ptr_with_raw_ptr_dispatch<T> d)
{
  using std::distance;
  if (distance(f, l) == 0)
    return d;
  d.wSP_->fence();
  using std::uninitialized_copy;
  if (d.wSP_->get_group().root())
    uninitialized_copy(f, l, to_address(d));
  d.wSP_->fence();
  mpi3::communicator(d.wSP_->get_group(), 0).barrier();
  return d + distance(f, l);
}

template<class T, class Size>
shm_ptr_with_raw_ptr_dispatch<T> uninitialized_default_construct_n(shm_ptr_with_raw_ptr_dispatch<T> f, Size n)
{
  if (n == 0)
    return f;
  f.wSP_->fence();
  if (f.wSP_->get_group().root())
  {
    T* current(to_address(f));
    try
    {
      for (; n > 0; ++current, --n)
        (::new ((void*)current) T());
    }
    catch (...)
    {
      throw;
    } // leak!
  }
  f.wSP_->fence();
  mpi3::communicator(f.wSP_->get_group(), 0).barrier();
  return f + n;
}

template<class T, class Size>
shm_ptr_with_raw_ptr_dispatch<T> uninitialized_value_construct_n(shm_ptr_with_raw_ptr_dispatch<T> f, Size n)
{
  if (n == 0)
    return f;
  f.wSP_->fence();
  if (f.wSP_->get_group().root())
  {
    T* current(to_address(f));
    try
    {
      for (; n > 0; ++current, --n)
        (::new ((void*)current) T());
    }
    catch (...)
    {
      throw;
    } // leak !!
  }
  f.wSP_->fence();
  mpi3::communicator(f.wSP_->get_group(), 0).barrier();
  return f + n;
}

template<class Alloc, class T, class Size>
shm_ptr_with_raw_ptr_dispatch<T> uninitialized_default_construct_n(Alloc& a, shm_ptr_with_raw_ptr_dispatch<T> f, Size n)
{
  return uninitialized_default_construct_n(f, n);
}

template<class Alloc, class T, class Size>
shm_ptr_with_raw_ptr_dispatch<T> uninitialized_value_construct_n(Alloc& a, shm_ptr_with_raw_ptr_dispatch<T> f, Size n)
{
  return uninitialized_value_construct_n(f, n);
}

template<class Alloc, class T, class Size>
shm_ptr_with_raw_ptr_dispatch<T> alloc_uninitialized_default_construct_n(Alloc& a,
                                                                         shm_ptr_with_raw_ptr_dispatch<T> f,
                                                                         Size n)
{
  return uninitialized_default_construct_n(f, n);
}

template<class Alloc, class T, class Size>
shm_ptr_with_raw_ptr_dispatch<T> alloc_uninitialized_value_construct_n(Alloc& a,
                                                                       shm_ptr_with_raw_ptr_dispatch<T> f,
                                                                       Size n)
{
  return uninitialized_value_construct_n(f, n);
}

} // namespace shm

namespace boost
{
namespace multi
{
template<typename T, typename Size, typename Q>
multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> fill_n(
    multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> first,
    Size n,
    Q const& val)
{
  if (n == 0)
    return first;
  base(first).wSP_->fence();
  if (first.wSP_->get_group().root())
  {
    auto current = first;
    try
    {
      for (; n > 0; ++current, --n)
        *(base(current)) = T(val);
    }
    catch (...)
    {
      throw;
    }
  }
  base(first).wSP_->fence();
  base(first).wSP_->fence();
  mpi3::communicator(base(first).wSP_->get_group(), 0).barrier();
  return first + n;
}

template<typename T, typename Q>
multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> fill_n(
    multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> first,
    multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> last,
    Q const& val)
{
  assert(stride(first) == stride(last));
  return fill_n(first, std::distance(first, last), val);
}


template<class Alloc, typename T, typename Size>
multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> uninitialized_fill_n(
    Alloc& a,
    multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> first,
    Size n,
    T const& val)
{
  if (n == 0)
    return first;
  base(first).wSP_->fence();
  if (first.wSP_->get_group().root())
  {
    auto current = first;
    try
    {
      for (; n > 0; ++current, --n)
        (::new ((void*)to_address(base(current))) T(val));
    }
    catch (...)
    {
      throw;
    }
  }
  base(first).wSP_->fence();
  base(first).wSP_->fence();
  mpi3::communicator(base(first).wSP_->get_group(), 0).barrier();
  return first + n;
}

template<class Alloc, typename T, typename Size>
multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> uninitialized_fill(
    Alloc& a,
    multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> first,
    multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> last,
    T const& val)
{
  assert(stride(first) == stride(last));
  return uninitialized_fill_n(a, first, std::distance(first, last), val);
}

template<class Alloc, typename T, typename Size>
multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> alloc_uninitialized_fill_n(
    Alloc& a,
    multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> first,
    Size n,
    T const& val)
{
  return uninitialized_fill_n(a, first, n, val);
}

template<class Alloc, typename T, typename Size>
multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> alloc_uninitialized_fill(
    Alloc& a,
    multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> first,
    multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> last,
    T const& val)
{
  assert(stride(first) == stride(last));
  return uninitialized_fill_n(a, first, std::distance(first, last), val);
}

template<class T, class Q1, class Q2, typename Size>
multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> copy_n(
    multi::array_iterator<Q1, 1, shm::shm_ptr_with_raw_ptr_dispatch<Q2>> first,
    Size n,
    multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> dest)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  if (n == 0)
    return dest;
  base(dest).wSP_->fence();
  base(first).wSP_->fence();
  if (base(dest).wSP_->get_group().root())
  {
    auto f = first;
    auto d = dest;
    for (; n > 0; ++f, ++d, --n)
      *d = *f;
  }
  base(dest).wSP_->fence();
  base(dest).wSP_->fence();
  mpi3::communicator(base(dest).wSP_->get_group(), 0).barrier();
  return dest + n;
}

template<class T, class ForwardIt, typename Size>
multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> copy_n(
    ForwardIt first,
    Size n,
    multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> dest)
{
  if (n == 0)
    return dest;
  base(dest).wSP_->fence();
  if (base(dest).wSP_->get_group().root())
  {
    auto f = first;
    auto d = dest;
    for (; n > 0; ++f, ++d, --n)
      *d = *f;
  }
  base(dest).wSP_->fence();
  base(dest).wSP_->fence();
  mpi3::communicator(base(dest).wSP_->get_group(), 0).barrier();
  return dest + n;
}

template<class T, class Q1, class Q2, typename Size>
multi::array_iterator<T, 1, T*> copy_n(multi::array_iterator<Q1, 1, shm::shm_ptr_with_raw_ptr_dispatch<Q2>> first,
                                       Size n,
                                       multi::array_iterator<T, 1, T*> dest)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  if (n == 0)
    return dest;
  base(first).wSP_->fence();
  {
    auto f = first;
    auto d = dest;
    for (; n > 0; ++f, ++d, --n)
      *d = *f;
  }
  return dest + n;
}

template<class T, class Q1, class Q2>
multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> copy(
    multi::array_iterator<Q1, 1, shm::shm_ptr_with_raw_ptr_dispatch<Q2>> first,
    multi::array_iterator<Q1, 1, shm::shm_ptr_with_raw_ptr_dispatch<Q2>> last,
    multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> dest)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  assert(stride(first) == stride(last));
  return copy_n(first, std::distance(first, last), dest);
}

template<class T, class ForwardIt>
multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> copy(
    ForwardIt first,
    ForwardIt last,
    multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> dest)
{
  assert(stride(first) == stride(last));
  return copy_n(first, std::distance(first, last), dest);
}

template<class T, class Q1, class Q2>
multi::array_iterator<T, 1, T*> copy(multi::array_iterator<Q1, 1, shm::shm_ptr_with_raw_ptr_dispatch<Q2>> first,
                                     multi::array_iterator<Q1, 1, shm::shm_ptr_with_raw_ptr_dispatch<Q2>> last,
                                     multi::array_iterator<T, 1, T*> dest)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  assert(stride(first) == stride(last));
  return copy_n(first, std::distance(first, last), dest);
}

template<class Alloc, class T, class Q, typename Size>
multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> uninitialized_copy_n(
    Alloc& a,
    multi::array_iterator<Q, 1, shm::shm_ptr_with_raw_ptr_dispatch<Q>> first,
    Size n,
    multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> dest)
{
  static_assert(std::is_same<typename std::decay<Q>::type, T>::value, "Wrong dispatch.\n");
  if (n == 0)
    return dest;
  base(first).wSP_->fence();
  base(dest).wSP_->fence();
  if (base(dest).wSP_->get_group().root())
  {
    auto f = first;
    auto d = dest;
    try
    {
      for (; n > 0; ++f, ++d, --n)
        (::new ((void*)to_address(base(d))) T(*f));
    }
    catch (...)
    {
      throw;
    }
  }
  base(dest).wSP_->fence();
  base(dest).wSP_->fence();
  mpi3::communicator(base(dest).wSP_->get_group(), 0).barrier();
  return dest + n;
}

template<class Alloc, class T, class Q, typename Size>
multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> alloc_uninitialized_copy_n(
    Alloc& a,
    multi::array_iterator<Q, 1, shm::shm_ptr_with_raw_ptr_dispatch<Q>> first,
    Size n,
    multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> dest)
{
  return uninitialized_copy_n(a, first, n, dest);
}


template<class Alloc, class T, class ForwardIt, typename Size>
multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> uninitialized_copy_n(
    Alloc& a,
    ForwardIt first,
    Size n,
    multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> dest)
{
  return copy_n(first, n, dest);
}

template<class Alloc, class T, class ForwardIt, typename Size>
multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> alloc_uninitialized_copy_n(
    Alloc& a,
    ForwardIt first,
    Size n,
    multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> dest)
{
  return copy_n(first, n, dest);
}

template<class Alloc, class T, class ForwardIt>
multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> uninitialized_copy(
    Alloc& a,
    ForwardIt first,
    ForwardIt last,
    multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> dest)
{
  assert(stride(first) == stride(last));
  return uninitialized_copy_n(a, first, std::distance(first, last), dest);
}

template<class Alloc, class T, class Q1, class Q2>
multi::array_iterator<T, 1, T*> uninitialized_copy(
    Alloc& a,
    multi::array_iterator<Q1, 1, shm::shm_ptr_with_raw_ptr_dispatch<Q2>> first,
    multi::array_iterator<Q1, 1, shm::shm_ptr_with_raw_ptr_dispatch<Q2>> last,
    multi::array_iterator<T, 1, T*> dest)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  if (std::distance(first, last) == 0)
    return dest;
  assert(stride(first) == stride(last));
  base(first).wSP_->fence();
  {
    auto d = dest;
    try
    {
      for (; first != last; ++first, ++d)
        (::new ((void*)to_address(base(d))) T(*first));
    }
    catch (...)
    {
      throw;
    }
  }
  mpi3::communicator(base(first).wSP_->get_group(), 0).barrier();
  return dest + std::distance(first, last);
}

template<class Alloc, class T, class ForwardIt>
multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> alloc_uninitialized_copy(
    Alloc& a,
    ForwardIt first,
    ForwardIt last,
    multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> dest)
{
  assert(stride(first) == stride(last));
  return uninitialized_copy_n(a, first, std::distance(first, last), dest);
}

template<class Alloc, class T, class Q1, class Q2>
multi::array_iterator<T, 1, T*> alloc_uninitialized_copy(
    Alloc& a,
    multi::array_iterator<Q1, 1, shm::shm_ptr_with_raw_ptr_dispatch<Q2>> first,
    multi::array_iterator<Q1, 1, shm::shm_ptr_with_raw_ptr_dispatch<Q2>> last,
    multi::array_iterator<T, 1, T*> dest)
{
  return uninitialized_copy(a, first, last, dest);
}

template<class Alloc, class T, class Size>
multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> uninitialized_default_construct_n(
    Alloc& a,
    multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> f,
    Size n)
{
  if (n == 0)
    return f;
  base(f).wSP_->fence();
  if (base(f).wSP_->get_group().root())
  {
    auto current(f);
    try
    {
      for (; n > 0; ++current, --n)
        (::new ((void*)to_address(base(current))) T());
      return current;
    }
    catch (...)
    {
      throw;
    }
  }
  base(f).wSP_->fence();
  mpi3::communicator(f.wSP_->get_group(), 0).barrier();
  return f + n;
}

template<class Alloc, class T, class Size>
multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> uninitialized_value_construct_n(
    Alloc& a,
    multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> f,
    Size n)
{
  return uninitialized_default_construct_n(a, f, n);
}

template<class Alloc, class T, class Size>
multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> alloc_uninitialized_default_construct_n(
    Alloc& a,
    multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> f,
    Size n)
{
  return uninitialized_default_construct_n(a, f, n);
}

template<class Alloc, class T, class Size>
multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> alloc_uninitialized_value_construct_n(
    Alloc& a,
    multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> f,
    Size n)
{
  return uninitialized_default_construct_n(a, f, n);
}

} // namespace multi
} // namespace boost

#endif
