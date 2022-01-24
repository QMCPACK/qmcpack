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

#ifndef AFQMC_DEVICE_POINTERS_HPP
#define AFQMC_DEVICE_POINTERS_HPP

#include <functional>
#include "Configuration.h"
#include <cassert>

#include "multi/array_ref.hpp"
#include "mpi3/shared_window.hpp"
#include "AFQMC/Memory/arch.hpp"
#include "AFQMC/Memory/raw_pointers.hpp"
#include "AFQMC/Numerics/device_kernels.hpp"

namespace device
{
using qmcplusplus::afqmc::to_address;

template<class T>
struct device_allocator;
template<class T>
struct device_pointer;
struct memory_resource;

// no const for now
template<class T>
struct device_reference
{
public:
  using value_type       = T;
  using decay_value_type = typename std::decay<T>::type;
  using pointer          = device_pointer<T>;

  // must construct through a gpu_ptr for now, to keep some sort of control/safety
  device_reference(pointer const& gpu_ptr) : impl_(to_address(gpu_ptr)) {}
  device_reference(device_reference<T> const& gpu_ref) = default;

  // assignment
  device_reference& operator=(device_reference const& x)
  {
    arch::memcopy(impl_, x.impl_, sizeof(T));
    return *this;
  }

  device_reference& operator=(value_type const& x)
  {
    arch::memcopy(impl_, std::addressof(x), sizeof(T));
    return *this;
  }

  template<class>
  friend struct device_reference;
  //friend class device_reference<typename std::decay<T>::type>;
  operator device_reference<T const>() const& { return device_reference<T const>{impl_}; }

  operator decay_value_type() const& { return this->val(); }

  // try getting rid of this
  //  operator decay_value_type&() & = delete;

  pointer operator&() const { return pointer{impl_}; }

  value_type* addressof() { return impl_; }
  value_type const* addressof() const { return impl_; }
  friend auto addressof(device_reference const& rhs) { return rhs.addressof(); }
  friend auto addressof(device_reference&& rhs) { return rhs.addressof(); }
  friend auto addressof(device_reference& rhs) { return rhs.addressof(); }

  void swap(device_reference& other) { std::swap(impl_, other.impl_); }

  device_reference& operator++(void)
  {
    kernels::op_plus(impl_, T(1));
    return *this;
  }
  device_reference& operator--(void)
  {
    kernels::op_minus(impl_, T(1));
    return *this;
  }

  value_type operator++(int)
  {
    value_type res;
    arch::memcopy(std::addressof(res), impl_, sizeof(T));
    return res++;
  }
  value_type operator--(int)
  {
    value_type res;
    arch::memcopy(std::addressof(res), impl_, sizeof(T));
    return res--;
  }

  device_reference& operator+=(const value_type& rhs)
  {
    kernels::op_plus(impl_, rhs);
    return *this;
  }
  device_reference& operator-=(const value_type& rhs)
  {
    kernels::op_minus(impl_, rhs);
    return *this;
  }
  device_reference& operator*=(const value_type& rhs)
  {
    kernels::op_times(impl_, rhs);
    return *this;
  }
  device_reference& operator/=(const value_type& rhs)
  {
    kernels::op_div(impl_, rhs);
    return *this;
  }

  device_reference& operator+=(device_reference const& rhs)
  {
    kernels::op_plus(impl_, rhs.val());
    return *this;
  }
  device_reference& operator-=(device_reference const& rhs)
  {
    kernels::op_minus(impl_, rhs.val());
    return *this;
  }
  device_reference& operator*=(device_reference const& rhs)
  {
    kernels::op_times(impl_, rhs.val());
    return *this;
  }
  device_reference& operator/=(device_reference const& rhs)
  {
    kernels::op_div(impl_, rhs.val());
    return *this;
  }

  friend value_type& operator+=(value_type& lhs, device_reference const& rhs)
  {
    lhs += rhs.val();
    return lhs;
  }
  friend value_type& operator-=(value_type& lhs, device_reference const& rhs)
  {
    lhs -= rhs.val();
    return lhs;
  }
  friend value_type& operator*=(value_type& lhs, device_reference const& rhs)
  {
    lhs *= rhs.val();
    return lhs;
  }
  friend value_type& operator/=(value_type& lhs, device_reference const& rhs)
  {
    lhs /= rhs.val();
    return lhs;
  }

  value_type operator+(value_type const& rhs) const { return this->val() + rhs; }
  value_type operator-(value_type const& rhs) const { return this->val() - rhs; }
  value_type operator/(value_type const& rhs) const { return this->val() / rhs; }
  value_type operator*(value_type const& rhs) const { return this->val() * rhs; }

  value_type operator+(device_reference const& rhs) const { return this->val() + rhs.val(); }
  value_type operator-(device_reference const& rhs) const { return this->val() - rhs.val(); }
  value_type operator*(device_reference const& rhs) const { return this->val() * rhs.val(); }
  value_type operator/(device_reference const& rhs) const { return this->val() / rhs.val(); }

  friend value_type operator+(value_type lhs, device_reference const& rhs) { return lhs + rhs.val(); }
  friend value_type operator-(value_type lhs, device_reference const& rhs) { return lhs - rhs.val(); }
  friend value_type operator*(value_type lhs, device_reference const& rhs) { return lhs * rhs.val(); }
  friend value_type operator/(value_type lhs, device_reference const& rhs) { return lhs / rhs.val(); }

  bool operator==(value_type const& rhs) const { return this->val() == rhs; }
  bool operator!=(value_type const& rhs) const { return this->val() != rhs; }
  bool operator>(value_type const& rhs) const { return this->val() > rhs; }
  bool operator<(value_type const& rhs) const { return this->val() < rhs; }
  bool operator>=(value_type const& rhs) const { return this->val() >= rhs; }
  bool operator<=(value_type const& rhs) const { return this->val() <= rhs; }

  bool operator==(device_reference const& rhs) const { return this->val() == rhs.val(); }
  bool operator!=(device_reference const& rhs) const { return this->val() != rhs.val(); }
  bool operator>(device_reference const& rhs) const { return this->val() > rhs.val(); }
  bool operator<(device_reference const& rhs) const { return this->val() < rhs.val(); }
  bool operator>=(device_reference const& rhs) const { return this->val() >= rhs.val(); }
  bool operator<=(device_reference const& rhs) const { return this->val() <= rhs.val(); }

  friend bool operator==(value_type const& lhs, device_reference const& rhs) { return lhs == rhs.val(); }
  friend bool operator!=(value_type const& lhs, device_reference const& rhs) { return lhs != rhs.val(); }
  friend bool operator>(value_type const& lhs, device_reference const& rhs) { return lhs > rhs.val(); }
  friend bool operator<(value_type const& lhs, device_reference const& rhs) { return lhs < rhs.val(); }
  friend bool operator>=(value_type const& lhs, device_reference const& rhs) { return lhs >= rhs.val(); }
  friend bool operator<=(value_type const& lhs, device_reference const& rhs) { return lhs <= rhs.val(); }

  friend std::ostream& operator<<(std::ostream& os, device_reference const& obj)
  {
    os << obj.val();
    return os;
  }

  friend std::istream& operator<<(std::istream& is, device_reference& obj)
  {
    value_type val;
    is >> val;
    obj = val;
    return is;
  }

  /*
  device_reference&  operator<<= (const T &rhs)
 
  device_reference&  operator>>= (const T &rhs)
 
  device_reference&  operator&= (const T &rhs)
 
  device_reference&  operator|= (const T &rhs)
 
  device_reference&  operator^= (const T &rhs)
*/

private:
  device_reference(T* ptr) : impl_(ptr) {}
  value_type* impl_;

  value_type val() const
  {
    decay_value_type res;
    arch::memcopy(std::addressof(res), impl_, sizeof(value_type));
    return value_type(res);
  }
};

struct base_device_pointer
{
  static arch::device_handles handles;
};

template<>
struct device_pointer<const void> : base_device_pointer
{
  using T            = const void;
  using element_type = void;
  using value_type   = void;
  const void* impl_;
  device_pointer(std::nullptr_t = nullptr) : impl_(nullptr) {}
  device_pointer(device_pointer const& other)            = default;
  device_pointer& operator=(device_pointer const& other) = default;
  device_pointer& operator=(std::nullptr_t)
  {
    impl_ = nullptr;
    return *this;
  }
  bool operator==(std::nullptr_t) const { return (impl_ == nullptr); }
  bool operator!=(std::nullptr_t) const { return not operator==(nullptr); }

private:
  device_pointer(T* impl__) : impl_(impl__) {}
  template<class>
  friend struct device_pointer;
};

template<>
struct device_pointer<void> : base_device_pointer
{
  using T            = void;
  using element_type = T;
  using value_type   = void;
  void* impl_;
  friend struct memory_resource;
  device_pointer(std::nullptr_t = nullptr) : impl_(nullptr) {}
  device_pointer(device_pointer const& other) = default;
  template<typename Q>
  device_pointer(device_pointer<Q> const& ptr) : impl_(reinterpret_cast<T*>(ptr.impl_))
  {}
  device_pointer& operator=(device_pointer const& other) = default;
  device_pointer& operator=(std::nullptr_t)
  {
    impl_ = nullptr;
    return *this;
  }
  bool operator==(std::nullptr_t) const { return (impl_ == nullptr); }
  bool operator!=(std::nullptr_t) const { return not operator==(nullptr); }
  explicit operator bool() const { return (impl_ != nullptr); }
  bool operator==(device_pointer const& other) const { return impl_ == other.impl_; }
  bool operator!=(device_pointer const& other) const { return not(*this == other); }
  bool operator<=(device_pointer<T> const& other) const { return impl_ <= other.impl_; }

private:
  device_pointer(T* impl__) : impl_(impl__) {}
  template<class>
  friend struct device_pointer;
};

// this class is not safe, since it allows construction of a gpu_ptr from a raw ptr
// which might not be in gpu memory. Fix this!!!
template<class T>
struct device_pointer : base_device_pointer
{
  using difference_type   = std::ptrdiff_t;
  using value_type        = std::decay_t<T>;
  using element_type      = T;
  using const_value_type  = T const;
  using pointer           = T*;
  using const_pointer     = T const*;
  using reference         = device_reference<T>;
  using const_reference   = device_reference<const T>;
  using iterator_category = std::random_access_iterator_tag;
  using rebind_const      = device_pointer<const T>;
  //  static const int memory_type = GPU_MEMORY_POINTER_TYPE;
  using default_allocator_type = device_allocator<value_type>;
  friend struct device_allocator<value_type>;
  friend struct device_pointer<typename std::decay<T>::type>;
  default_allocator_type default_allocator() const { return device_allocator<value_type>{}; };
  device_pointer() = default;
  device_pointer(std::nullptr_t) : impl_(nullptr) {}
  template<typename Q>
  device_pointer(device_pointer<Q> const& ptr) : impl_(reinterpret_cast<T*>(ptr.impl_))
  {}
  reference operator*() const { return reference(device_pointer{impl_}); }
  reference operator[](std::ptrdiff_t n) const { return reference(device_pointer{impl_ + n}); }
  // is this correct? should I just delete it?
  T* operator->() const { return impl_; }
  explicit operator bool() const { return (impl_ != nullptr); }
  auto operator+(std::ptrdiff_t n) const { return device_pointer{impl_ + n}; }
  auto operator-(std::ptrdiff_t n) const { return device_pointer{impl_ - n}; }
  std::ptrdiff_t operator-(device_pointer other) const { return std::ptrdiff_t(impl_ - other.impl_); }
  operator device_pointer<T const>() const { return device_pointer<T const>{impl_}; }
  operator device_pointer<void const>() const { return device_pointer<void const>{impl_}; }
  device_pointer& operator++()
  {
    ++impl_;
    return *this;
  }
  device_pointer& operator--()
  {
    --impl_;
    return *this;
  }
  device_pointer& operator+=(std::ptrdiff_t d)
  {
    impl_ += d;
    return *this;
  }
  device_pointer& operator-=(std::ptrdiff_t d)
  {
    impl_ -= d;
    return *this;
  }
  bool operator==(device_pointer const& other) const { return impl_ == other.impl_; }
  bool operator!=(device_pointer const& other) const { return not(*this == other); }
  bool operator<=(device_pointer<T> const& other) const { return impl_ <= other.impl_; }

  T* to_address() const { return impl_; }
  friend decltype(auto) to_address(device_pointer const& self) { return self.to_address(); }
  template<class Q>
  device_pointer<Q> pointer_cast()
  {
    device_pointer<Q> res;
    res.impl_ = reinterpret_cast<Q*>(impl_);
    return res;
  }
  template<class Q>
  friend device_pointer<Q> pointer_cast(device_pointer&& self)
  {
    return self.pointer_cast<Q>();
    //    return device_pointer<Q>{reinterpret_cast<Q*>(self.impl_)};
  }
  T* impl_;

private:
  device_pointer(T* impl__) : impl_(impl__) {}
};

//static size_t TotalGPUAlloc=0;

/*
 * Incomplete!!! Need to fix construct and destroy
 */
template<class T>
struct device_allocator
{
  template<class U>
  struct rebind
  {
    using other = device_allocator<U>;
  };
  using element_type     = T;
  using value_type       = T;
  using const_value_type = T const;
  using pointer          = device_pointer<T>;
  using const_pointer    = device_pointer<T const>;
  using reference        = device_reference<T>;
  using const_reference  = device_reference<T const>;
  using size_type        = std::size_t;
  using difference_type  = std::ptrdiff_t;

  device_allocator()                              = default;
  ~device_allocator()                             = default;
  device_allocator(device_allocator const& other) = default;
  template<class U>
  device_allocator(device_allocator<U> const& other)
  {}

  device_pointer<T> allocate(size_type n, const void* hint = 0)
  {
    if (n == 0)
      return device_pointer<T>{};
    T* p;
    arch::malloc((void**)&p, n * sizeof(T));
    device_pointer<T> res;
    res.impl_ = p;
    return res;
  }
  void deallocate(device_pointer<T> ptr, size_type) { arch::free(ptr.impl_); }
  bool operator==(device_allocator const& other) const { return true; }
  bool operator!=(device_allocator const& other) const { return false; }
  template<class U, class... Args>
  void construct(U p, Args&&... args)
  {
    static_assert(std::is_trivially_copy_constructible<value_type>{},
                  "!"); // ::new((void*)p) U(std::forward<Args>(args)...);
  }
  template<class U>
  void destroy(U p)
  {
    static_assert(std::is_trivially_destructible<value_type>{}, "!"); // p->~U();
  }
  template<class InputIt, class ForwardIt>
  ForwardIt alloc_uninitialized_copy(InputIt first, InputIt last, ForwardIt d_first)
  {
    static_assert(std::is_trivially_copy_constructible<value_type>{}, "!");
    static_assert(std::is_trivially_destructible<value_type>{}, "!");
    std::advance(d_first, std::distance(first, last));
    return d_first;
  }
};

struct memory_resource
{
  using pointer = device_pointer<void>;

  device_pointer<void> allocate(std::size_t size, std::size_t alignment = alignof(std::max_align_t))
  {
    if (size == 0)
      return device_pointer<void>{};
    void* p;
    arch::malloc(&p, size);
    return device_pointer<void>{p};
  }

  void deallocate(device_pointer<void> ptr, std::size_t = alignof(std::max_align_t)) { arch::free(ptr.impl_); }

  bool operator==(memory_resource const& other) const { return true; }

  bool operator!=(memory_resource const& other) const { return false; }
};

// finish later
template<class T>
struct constructor
{
  template<class U>
  struct rebind
  {
    using other = constructor<U>;
  };
  using value_type         = T;
  using pointer            = device_pointer<T>;
  using const_pointer      = device_pointer<T const>;
  using void_pointer       = device_pointer<void>;
  using const_void_pointer = device_pointer<void const>;
  using size_type          = std::size_t;
  using difference_type    = std::ptrdiff_t;

  constructor()                         = default;
  ~constructor()                        = default;
  constructor(constructor const& other) = default;
  template<class U>
  constructor(constructor<U> const& other)
  {}

  template<class U, class... Args>
  void construct(U p, Args&&... args)
  {}

  template<class U>
  void destroy(U p)
  {}
};

/* Don't know how to implement this on the kernel side, without propagating the 
 * cuda code upstream due to the template needed to pass the UnaryOperator
template<class T, class F>
F for_each(device_pointer<T> first, device_pointer<T> last, F f){
        if(first == last) return f;
        return kernels::for_each(to_address(first), to_address(last), f);
}
*/
/**************** copy_n *****************/
/*
template<typename T, typename Size>
device_pointer<T> copy_n(device_pointer<T const> const A, Size n, device_pointer<T> B) {
  arch::memcopy(to_address(B),to_address(A),n*sizeof(T)); 
  return B+n;
}
*/

template<typename T, typename ForwardIt, typename Size>
device_pointer<T> copy_n(ForwardIt A, Size n, device_pointer<T> B)
{
  arch::memcopy(to_address(B), to_address(A), n * sizeof(T));
  return B + n;
}

/*
// MAM: removing this to force copy to other pointers from the other ptr point of view
template<typename T, typename ForwardIt, typename Size>
ForwardIt copy_n(device_pointer<T const> const A, Size n, ForwardIt B) {
  arch::memcopy(to_address(B),to_address(A),n*sizeof(T));
  return B+n;
}
*/

template<typename Q, typename T, typename Size>
device_pointer<T> copy_n(Q* const A, Size n, device_pointer<T> B)
{
  static_assert(std::is_same<typename std::decay<Q>::type, T>::value, "Wrong dispatch.\n");
  arch::memcopy(to_address(B), A, n * sizeof(T));
  return B + n;
}

template<typename T, typename Q, typename Size>
T* copy_n(device_pointer<Q> const A, Size n, T* B)
{
  static_assert(std::is_same<typename std::decay<Q>::type, T>::value, "Wrong dispatch.\n");
  arch::memcopy(B, to_address(A), n * sizeof(T));
  return B + n;
}

/**************** copy *****************/
/*
template<typename T>
device_pointer<T> copy(device_pointer<T const> const Abeg, 
                       device_pointer<T const> const Aend, 
                       device_pointer<T> B) {
  return copy_n(Abeg,std::distance(Abeg,Aend),B);
}
*/

template<typename T, typename ForwardIt>
device_pointer<T> copy(ForwardIt const Abeg, ForwardIt const Aend, device_pointer<T> B)
{
  return copy_n(Abeg, std::distance(Abeg, Aend), B);
}

template<typename T, typename Q>
device_pointer<T> copy(Q* const Abeg, Q* const Aend, device_pointer<T> B)
{
  static_assert(std::is_same<typename std::decay<Q>::type, T>::value, "Wrong dispatch.\n");
  return copy_n(Abeg, std::distance(Abeg, Aend), B);
}

template<typename T, typename Q>
T* copy(device_pointer<Q> const Abeg, device_pointer<Q> const Aend, T* B)
{
  static_assert(std::is_same<typename std::decay<Q>::type, T>::value, "Wrong dispatch.\n");
  return copy_n(Abeg, std::distance(Abeg, Aend), B);
}

/*
template<typename T, typename ForwardIt>
ForwardIt copy(device_pointer<T const> const Abeg, device_pointer<T const> const Aend, ForwardIt B) {
  copy_n(Abeg,std::distance(Abeg,Aend),to_address(B));
  return B + std::distance(Abeg,Aend);  
}
*/

/**************** copy_n_cast *****************/
// NOTE: Eliminate this routine, merge with copy_n and dispatch to kernel call
// if types of pointers are not the same (without cv qualifiers)!!!
template<typename T, typename Q, typename Size>
device_pointer<Q> copy_n_cast(device_pointer<T> const A, Size n, device_pointer<Q> B)
{
  kernels::copy_n_cast(to_address(A), n, to_address(B));
  return B + n;
}

template<typename T, typename Size>
device_pointer<T> copy_n_cast(device_pointer<T> const A, Size n, device_pointer<T> B)
{
  return copy_n(A, n, B);
}

template<typename T, typename Q, typename Size>
device_pointer<Q> copy_n_cast(T* const A, Size n, device_pointer<Q> B)
{
  throw std::runtime_error(" Error: copy_n_cast(T*,n,gpu_ptr) is disabled.");
  return B + n;
}

template<typename T, typename Q, typename Size>
Q* copy_n_cast(device_pointer<T> const A, Size n, Q* B)
{
  throw std::runtime_error(" Error: copy_n_cast(gpu_ptr,n,T*) is disabled.");
  return B + n;
}

/**************** inplace_cast *****************/
template<typename T, typename Q, typename Size>
void inplace_cast(device_pointer<T> A, Size n)
{
  T* A_(to_address(A));
  Q* B_(reinterpret_cast<Q*>(A_));
  kernels::inplace_cast(n, A_, B_);
}

/**************** fill_n *****************/
//template<typename T, typename Size, typename... Args>
//device_pointer<T> fill_n(device_pointer<T> first, Size n, Args&&...args){
//  if(n == 0) return first;
//  kernels::fill_n(to_address(first), n, std::forward<Args>(args)...);
//  return first + n;
//}

template<typename T, typename Size, typename Q>
device_pointer<T> fill_n(device_pointer<T> first, Size n, Q const& val)
{
  if (n == 0)
    return first;
  kernels::fill_n(to_address(first), n, T(val));
  return first + n;
}

template<typename T, typename Q>
device_pointer<T> fill(device_pointer<T> first, device_pointer<T> last, Q const& val)
{
  return fill_n(first, std::distance(first, last), T(val));
}


/**************** uninitialized_fill_n *****************/

template<typename T, typename Size, typename Q>
device_pointer<T> uninitialized_fill_n(device_pointer<T> first, Size n, Q const& val)
{
  if (n == 0)
    return first;
  //kernels::uninitialized_fill_n(to_address(first), n, val);
  kernels::fill_n(to_address(first), n, T(val));
  return first + n;
}

template<typename T, typename Q>
device_pointer<T> uninitialized_fill(device_pointer<T> first, device_pointer<T> last, Q const& val)
{
  return uninitialized_fill_n(first, std::distance(first, last), T(val));
}

template<class Alloc, typename T, typename Size, typename Q>
device_pointer<T> uninitialized_fill_n(Alloc& a, device_pointer<T> first, Size n, Q const& val)
{
  if (n == 0)
    return first;
  //kernels::uninitialized_fill_n(to_address(first), n, val);
  kernels::fill_n(to_address(first), n, T(val));
  return first + n;
}

template<class Alloc, typename T, typename Q>
device_pointer<T> uninitialized_fill(Alloc& a, device_pointer<T> first, device_pointer<T> last, Q const& val)
{
  return uninitialized_fill_n(a, first, std::distance(first, last), T(val));
}

template<class Alloc, typename T, typename Size, typename Q>
device_pointer<T> alloc_uninitialized_fill_n(Alloc& a, device_pointer<T> first, Size n, Q const& val)
{
  if (n == 0)
    return first;
  //kernels::uninitialized_fill_n(to_address(first), n, val);
  kernels::fill_n(to_address(first), n, T(val));
  return first + n;
}

template<class Alloc, typename T, typename Q>
device_pointer<T> alloc_uninitialized_fill(Alloc& a, device_pointer<T> first, device_pointer<T> last, Q const& val)
{
  return uninitialized_fill_n(a, first, std::distance(first, last), T(val));
}

/******************/

template<typename T, typename Size>
device_pointer<T> uninitialized_default_construct_n(device_pointer<T> first, Size n)
{
  return uninitialized_fill_n(first, n, T());
}

template<typename T>
device_pointer<T> uninitialized_default_construct(device_pointer<T> first, device_pointer<T> last)
{
  return uninitialized_fill_n(first, std::distance(first, last), T());
}

template<class Alloc, typename T, typename Size>
device_pointer<T> uninitialized_default_construct_n(Alloc& a, device_pointer<T> first, Size n)
{
  return uninitialized_fill_n(first, n, T());
}

template<class Alloc, typename T>
device_pointer<T> uninitialized_default_construct(Alloc& a, device_pointer<T> first, device_pointer<T> last)
{
  return uninitialized_fill_n(first, std::distance(first, last), T());
}

template<class Alloc, typename T, typename Size>
device_pointer<T> alloc_uninitialized_default_construct_n(Alloc& a, device_pointer<T> first, Size n)
{
  return uninitialized_fill_n(first, n, T());
}

template<class Alloc, typename T>
device_pointer<T> alloc_uninitialized_default_construct(Alloc& a, device_pointer<T> first, device_pointer<T> last)
{
  return uninitialized_fill_n(first, std::distance(first, last), T());
}

template<typename T, typename Size>
device_pointer<T> uninitialized_value_construct_n(device_pointer<T> first, Size n)
{
  return uninitialized_fill_n(first, n, T());
}

template<typename T>
device_pointer<T> uninitialized_value_construct_n(device_pointer<T> first, device_pointer<T> last)
{
  return uninitialized_fill_n(first, std::distance(first, last), T());
}

template<class Alloc, typename T, typename Size>
device_pointer<T> uninitialized_value_construct_n(Alloc& a, device_pointer<T> first, Size n)
{
  return uninitialized_fill_n(first, n, T());
}

template<class Alloc, typename T>
device_pointer<T> uninitialized_value_construct_n(Alloc& a, device_pointer<T> first, device_pointer<T> last)
{
  return uninitialized_fill_n(first, std::distance(first, last), T());
}

template<class Alloc, typename T, typename Size>
device_pointer<T> alloc_uninitialized_value_construct_n(Alloc& a, device_pointer<T> first, Size n)
{
  return uninitialized_fill_n(first, n, T());
}

template<class Alloc, typename T>
device_pointer<T> alloc_uninitialized_value_construct_n(Alloc& a, device_pointer<T> first, device_pointer<T> last)
{
  return uninitialized_fill_n(first, std::distance(first, last), T());
}

/**************** uninitialized_copy_n *****************/
/*
template<typename T, typename Size> 
device_pointer<T> uninitialized_copy_n(device_pointer<T> first, Size n, device_pointer<T> dest){
  if(n == 0) return dest;
  kernels::uninitialized_copy_n(n,to_address(first), 1, to_address(dest), 1);
  return dest + n;
}

template<class T> 
device_pointer<T> uninitialized_copy(device_pointer<T> first, device_pointer<T> last, device_pointer<T> dest){
  return uninitialized_copy_n(first,std::distance(first,last),dest); 
}
*/
// only trivial types for now, no placement new yet
template<typename T, typename Size>
device_pointer<T> uninitialized_copy_n(device_pointer<T> A, Size n, device_pointer<T> B)
{
  return copy_n(A, n, B);
}

template<typename T, typename Size, typename Q>
device_pointer<T> uninitialized_copy_n(T* A, Size n, device_pointer<Q> B)
{
  return copy_n(A, n, B);
}

template<typename T, typename Size, typename Q>
T* uninitialized_copy_n(device_pointer<T> A, Size n, Q* B)
{
  static_assert(std::is_trivially_assignable<Q&, T>{}, "!");
  return copy_n(A, n, B);
}

template<class Alloc, typename T, typename Size>
device_pointer<T> uninitialized_copy_n(Alloc& a, device_pointer<T> A, Size n, device_pointer<T> B)
{
  return copy_n(A, n, B);
}

template<class Alloc, typename T, typename Size, typename Q>
device_pointer<T> alloc_uninitialized_copy_n(Alloc& a, T* A, Size n, device_pointer<Q> B)
{
  return copy_n(A, n, B);
}

template<class Alloc, typename T, typename Size, typename Q>
T* alloc_uninitialized_copy_n(Alloc& a, device_pointer<T> A, Size n, Q* B)
{
  return copy_n(A, n, B);
}

template<class Alloc, typename T, typename Size>
device_pointer<T> alloc_uninitialized_copy_n(Alloc& a, device_pointer<T> A, Size n, device_pointer<T> B)
{
  return copy_n(A, n, B);
}

template<class Alloc, typename T, typename Size>
device_pointer<T> alloc_uninitialized_copy_n(Alloc& a, T const* A, Size n, device_pointer<T> B)
{
  return copy_n(A, n, B);
}

template<class Alloc, class Q, typename T, typename Size>
device_pointer<T> alloc_uninitialized_copy_n(Alloc& a, Q A, Size n, device_pointer<T> B)
{
  return copy_n(A, n, B);
}

template<class Alloc, typename T, typename Size>
T* alloc_uninitialized_copy_n(Alloc& a, device_pointer<T> A, Size n, T* B)
{
  return copy_n(A, n, B);
}

template<class Alloc, typename T, typename Size>
device_pointer<T> alloc_uninitialized_move_n(Alloc& a, device_pointer<T> A, Size n, device_pointer<T> B)
{
  return copy_n(A, n, B);
}

template<class Alloc, typename T, typename Size>
device_pointer<T> alloc_uninitialized_move_n(Alloc& a, T const* A, Size n, device_pointer<T> B)
{
  return copy_n(A, n, B);
}

template<class Alloc, class Q, typename T, typename Size>
device_pointer<T> alloc_uninitialized_move_n(Alloc& a, Q A, Size n, device_pointer<T> B)
{
  return copy_n(A, n, B);
}

template<class Alloc, typename T, typename Size>
T* alloc_uninitialized_move_n(Alloc& a, device_pointer<T> A, Size n, T* B)
{
  return copy_n(A, n, B);
}

/**************** uninitialized_copy *****************/
template<typename T>
device_pointer<T> uninitialized_copy(device_pointer<T> const Abeg, device_pointer<T> const Aend, device_pointer<T> B)
{
  return copy_n(Abeg, std::distance(Abeg, Aend), B);
}

template<typename T, typename Q>
device_pointer<T> uninitialized_copy(T* Abeg, T* Aend, device_pointer<Q> B)
{
  return copy_n(Abeg, std::distance(Abeg, Aend), B);
}

template<typename T>
T* uninitialized_copy(device_pointer<T> const Abeg, device_pointer<T> const Aend, T* B)
{
  return copy_n(Abeg, std::distance(Abeg, Aend), B);
}

template<class Alloc, typename T>
device_pointer<T> uninitialized_copy(Alloc& a,
                                     device_pointer<T> const Abeg,
                                     device_pointer<T> const Aend,
                                     device_pointer<T> B)
{
  return copy_n(Abeg, std::distance(Abeg, Aend), B);
}

template<class Alloc, typename T>
device_pointer<T> uninitialized_copy(Alloc& a, T* const Abeg, T* const Aend, device_pointer<T> B)
{
  return copy_n(Abeg, std::distance(Abeg, Aend), B);
}

template<class Alloc, typename T>
T* uninitialized_copy(Alloc& a, device_pointer<T> const Abeg, device_pointer<T> const Aend, T* B)
{
  return copy_n(Abeg, std::distance(Abeg, Aend), B);
}

template<class Alloc, typename T>
device_pointer<T> alloc_uninitialized_copy(Alloc& a,
                                           device_pointer<T> const Abeg,
                                           device_pointer<T> const Aend,
                                           device_pointer<T> B)
{
  return copy_n(Abeg, std::distance(Abeg, Aend), B);
}

template<class Alloc, typename T, typename Q>
device_pointer<Q> uninitialized_copy(T* Abeg, T* Aend, device_pointer<Q> B)
{
  static_assert(std::is_trivially_assignable<Q&, T>{}, "!");
  return copy_n(Abeg, std::distance(Abeg, Aend), B);
}

template<class Alloc, typename T, typename Q>
device_pointer<Q> alloc_uninitialized_copy(Alloc& a, T* Abeg, T* Aend, device_pointer<Q> B)
{
  static_assert(std::is_trivially_assignable<Q&, T>{}, "!");
  return uninitialized_copy(Abeg, Aend, B);
}

template<class Alloc, typename T>
T* alloc_uninitialized_copy(Alloc& a, device_pointer<T> const Abeg, device_pointer<T> const Aend, T* B)
{
  return copy_n(Abeg, std::distance(Abeg, Aend), B);
}

/**************** destroy_n *****************/
// NOTE: Not sure what to do here
// should at least guard against non-trivial types
template<typename T, typename Size>
device_pointer<T> destroy_n(device_pointer<T> first, Size n)
{
  return first + n;
}

template<class Alloc, typename T, typename Size>
device_pointer<T> destroy_n(Alloc& a, device_pointer<T> first, Size n)
{
  return first + n;
}

template<class Alloc, typename T, typename Size>
device_pointer<T> alloc_destroy_n(Alloc& a, device_pointer<T> first, Size n)
{
  return first + n;
}

/**************** print *****************/
template<typename T>
void print(std::string str, device_pointer<T> p, int n)
{
  kernels::print(str, to_address(p), n);
}

template<typename T, typename Q>
void fill2D(int n, int m, device::device_pointer<T> first, int lda, Q const& val)
{
  assert(lda >= m);
  kernels::fill2D_n(n, m, to_address(first), lda, T(val));
}

} // namespace device

namespace boost
{
namespace multi
{
//using qmcplusplus::afqmc::to_address;

/**************** copy *****************/
// Can always call cudaMemcopy2D like you do in the blas backend

template<typename T, typename Size>
multi::array_iterator<T, 1, device::device_pointer<T>> fill_n(
    multi::array_iterator<T, 1, device::device_pointer<T>> first,
    Size n,
    T const& val)
{
  if (n == 0)
    return first;
  kernels::fill_n(to_address(base(first)), n, stride(first), val);
  return first + n;
}

template<typename T, typename Size>
multi::array_iterator<T, 1, device::device_pointer<T>> fill(
    multi::array_iterator<T, 1, device::device_pointer<T>> first,
    multi::array_iterator<T, 1, device::device_pointer<T>> last,
    T const& val)
{
  assert(stride(first) == stride(last));
  if (std::distance(first, last) == 0)
    return first;
  kernels::fill_n(to_address(base(first)), std::distance(first, last), stride(first), val);
  return first + std::distance(first, last);
}

template<class Alloc, typename T, typename Size>
multi::array_iterator<T, 1, device::device_pointer<T>> uninitialized_fill_n(
    Alloc& a,
    multi::array_iterator<T, 1, device::device_pointer<T>> first,
    Size n,
    T const& val)
{
  if (n == 0)
    return first;
  kernels::fill_n(to_address(base(first)), n, stride(first), val);
  return first + n;
}

template<class Alloc, typename T, typename Size>
multi::array_iterator<T, 1, device::device_pointer<T>> uninitialized_fill(
    Alloc& a,
    multi::array_iterator<T, 1, device::device_pointer<T>> first,
    multi::array_iterator<T, 1, device::device_pointer<T>> last,
    T const& val)
{
  assert(stride(first) == stride(last));
  if (std::distance(first, last) == 0)
    return first;
  kernels::fill_n(to_address(base(first)), std::distance(first, last), stride(first), val);
  return first + std::distance(first, last);
}

template<class Alloc, typename T, typename Size>
multi::array_iterator<T, 1, device::device_pointer<T>> alloc_uninitialized_fill_n(
    Alloc& a,
    multi::array_iterator<T, 1, device::device_pointer<T>> first,
    Size n,
    T const& val)
{
  if (n == 0)
    return first;
  kernels::fill_n(to_address(base(first)), n, stride(first), val);
  return first + n;
}

template<class Alloc, typename T, typename Size>
multi::array_iterator<T, 1, device::device_pointer<T>> alloc_uninitialized_fill(
    Alloc& a,
    multi::array_iterator<T, 1, device::device_pointer<T>> first,
    multi::array_iterator<T, 1, device::device_pointer<T>> last,
    T const& val)
{
  assert(stride(first) == stride(last));
  if (std::distance(first, last) == 0)
    return first;
  kernels::fill_n(to_address(base(first)), std::distance(first, last), stride(first), val);
  return first + std::distance(first, last);
}

template<class T, class Q1, class Q2>
multi::array_iterator<T, 1, device::device_pointer<T>> copy(
    multi::array_iterator<Q1, 1, device::device_pointer<Q2>> first,
    multi::array_iterator<Q1, 1, device::device_pointer<Q2>> last,
    multi::array_iterator<T, 1, device::device_pointer<T>> dest)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  assert(stride(first) == stride(last));
  if (std::distance(first, last) == 0)
    return dest;
  arch::memcopy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                  sizeof(T), std::distance(first, last), arch::memcopyD2D);
  return dest + std::distance(first, last);
}

template<class T, class ForwardIt>
multi::array_iterator<T, 1, device::device_pointer<T>> copy(ForwardIt first,
                                                            ForwardIt last,
                                                            multi::array_iterator<T, 1, device::device_pointer<T>> dest)
{
  assert(stride(first) == stride(last));
  if (std::distance(first, last) == 0)
    return dest;
  using qmcplusplus::afqmc::to_address;
  arch::memcopy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                  sizeof(T), std::distance(first, last));
  return dest + std::distance(first, last);
}

template<typename T, typename Q, typename QQ>
multi::array_iterator<T, 1, device::device_pointer<T>> uninitialized_copy(
    T* first,
    T* last,
    multi::array_iterator<Q, 1, device::device_pointer<QQ>> dest)
{
  static_assert(std::is_trivially_assignable<QQ&, T>{}, "!");
  assert(stride(first) == stride(last));
  if (std::distance(first, last) == 0)
    return dest;
  using qmcplusplus::afqmc::to_address;
  arch::memcopy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * 1, sizeof(T),
                  std::distance(first, last));
  return dest + std::distance(first, last);
}

template<class ForwardIt, class Q1, class Q2>
ForwardIt copy(multi::array_iterator<Q1, 1, device::device_pointer<Q2>> first,
               multi::array_iterator<Q1, 1, device::device_pointer<Q2>> last,
               ForwardIt dest)
{
  using T = typename std::decay<typename ForwardIt::value_type>::type;
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  assert(stride(first) == stride(last));
  if (std::distance(first, last) == 0)
    return dest;
  using qmcplusplus::afqmc::to_address;
  arch::memcopy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                  sizeof(T), std::distance(first, last), arch::memcopyD2H);
  return dest + std::distance(first, last);
}

template<class T, class Q1, class Q2, typename Size>
multi::array_iterator<T, 1, device::device_pointer<T>> copy_n(
    multi::array_iterator<Q1, 1, device::device_pointer<Q2>> first,
    Size N,
    multi::array_iterator<T, 1, device::device_pointer<T>> dest)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  if (N == 0)
    return dest;
  arch::memcopy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                  sizeof(T), N, arch::memcopyD2D);
  return dest + N;
}

template<class T, class ForwardIt, typename Size>
multi::array_iterator<T, 1, device::device_pointer<T>> copy_n(
    ForwardIt first,
    Size n,
    multi::array_iterator<T, 1, device::device_pointer<T>> dest)
{
  if (n == 0)
    return dest;
  using qmcplusplus::afqmc::to_address;
  arch::memcopy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                  sizeof(T), n);
  return dest + n;
}

template<class ForwardIt, class Q1, class Q2, typename Size>
ForwardIt copy_n(multi::array_iterator<Q1, 1, device::device_pointer<Q2>> first, Size N, ForwardIt dest)
{
  using T = typename std::decay<typename ForwardIt::value_type>::type;
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  if (N == 0)
    return dest;
  using qmcplusplus::afqmc::to_address;
  arch::memcopy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                  sizeof(T), N, arch::memcopyD2H);
  return dest + N;
}

template<class Q, class QQ, class T, class TT>
multi::array_iterator<T, 1, device::device_pointer<T>> uninitialized_copy(
    multi::array_iterator<Q, 1, device::device_pointer<QQ>> first,
    multi::array_iterator<Q, 1, device::device_pointer<QQ>> last,
    multi::array_iterator<T, 1, device::device_pointer<TT>> dest)
{
  static_assert(std::is_trivially_assignable<TT&, QQ&>{}, "!");
  assert(stride(first) == stride(last));
  if (std::distance(first, last) == 0)
    return dest;
  using qmcplusplus::afqmc::to_address;
  arch::memcopy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                  sizeof(T), std::distance(first, last));
  return dest + std::distance(first, last);
}

template<class Alloc, class T, class ForwardIt>
multi::array_iterator<T, 1, device::device_pointer<T>> alloc_uninitialized_copy(
    Alloc& a,
    ForwardIt first,
    ForwardIt last,
    multi::array_iterator<T, 1, device::device_pointer<T>> dest)
{
  assert(stride(first) == stride(last));
  if (std::distance(first, last) == 0)
    return dest;
  using qmcplusplus::afqmc::to_address;
  arch::memcopy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                  sizeof(T), std::distance(first, last));
  return dest + std::distance(first, last);
}

template<class Alloc, class Q, class QQ, class T, class TT>
multi::array_iterator<T, 1, device::device_pointer<T>> alloc_uninitialized_copy(
    Alloc& a,
    multi::array_iterator<Q, 1, device::device_pointer<QQ>> first,
    multi::array_iterator<Q, 1, device::device_pointer<QQ>> last,
    multi::array_iterator<T, 1, device::device_pointer<TT>> dest)
{
  static_assert(std::is_trivially_assignable<TT&, QQ&>{}, "!");
  assert(stride(first) == stride(last));
  if (std::distance(first, last) == 0)
    return dest;
  using qmcplusplus::afqmc::to_address;
  arch::memcopy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                  sizeof(T), std::distance(first, last));
  return dest + std::distance(first, last);
}

/*
template<class Alloc, class T, class Q>
multi::array_iterator<T, 1, device::device_pointer<T>> uninitialized_copy( 
                         Alloc &a,
                         multi::array_iterator<Q, 1, device::device_pointer<Q>> first,
                         multi::array_iterator<Q, 1, device::device_pointer<Q>> last,
                         multi::array_iterator<T, 1, device::device_pointer<T>> dest ){
  static_assert(std::is_same<typename std::decay<Q>::type,T>::value,"Wrong dispatch.\n");
  assert( stride(first) == stride(last) );
  if(std::distance(first,last) == 0 ) return dest;
  if(cudaSuccess != cudaMemcpy2D(to_address(base(dest)),sizeof(T)*stride(dest),
                                 to_address(base(first)),sizeof(T)*stride(first),
                                 sizeof(T),std::distance(first,last),cudaMemcpyDeviceToDevice))
      throw std::runtime_error("Error: cudaMemcpy2D returned error code.");
  return dest+std::distance(first,last); 
}

template<class Alloc, class T, class Q>
multi::array_iterator<T, 1, device::device_pointer<T>> uninitialized_copy(
                         Alloc &a,
                         multi::array_iterator<Q, 1, boost::mpi3::intranode::array_ptr<Q>> first,
                         multi::array_iterator<Q, 1, boost::mpi3::intranode::array_ptr<Q>> last,
                         multi::array_iterator<T, 1, device::device_pointer<T>> dest ){
  static_assert(std::is_same<typename std::decay<Q>::type,T>::value,"Wrong dispatch.\n");
  assert( stride(first) == stride(last) );
  if(std::distance(first,last) == 0 ) return dest;
  if(cudaSuccess != cudaMemcpy2D(to_address(base(dest)),sizeof(T)*stride(dest),
                                 to_address(base(first)),sizeof(T)*stride(first),
                                 sizeof(T),std::distance(first,last),cudaMemcpyHostToDevice))
      throw std::runtime_error("Error: cudaMemcpy2D returned error code.");
  return dest+std::distance(first,last);
}

template<class Alloc, class T, class Q1, class Q2>
multi::array_iterator<T, 1, device::device_pointer<T>> uninitialized_copy(
                         Alloc &a,
                         multi::array_iterator<Q1, 1, Q2*> first,
                         multi::array_iterator<Q1, 1, Q2*> last,
                         multi::array_iterator<T, 1, device::device_pointer<T>> dest ){
  static_assert(std::is_same<typename std::decay<Q1>::type,T>::value,"Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type,T>::value,"Wrong dispatch.\n");
  assert( stride(first) == stride(last) );
  if(std::distance(first,last) == 0 ) return dest;
  if(cudaSuccess != cudaMemcpy2D(to_address(base(dest)),sizeof(T)*stride(dest),
                                 to_address(base(first)),sizeof(T)*stride(first),
                                 sizeof(T),std::distance(first,last),cudaMemcpyHostToDevice))
      throw std::runtime_error("Error: cudaMemcpy2D returned error code.");
  return dest+std::distance(first,last);
}
*/

template<class Alloc, class T, class Q1, class Q2>
multi::array_iterator<T, 1, T*> uninitialized_copy(Alloc& a,
                                                   multi::array_iterator<Q1, 1, device::device_pointer<Q2>> first,
                                                   multi::array_iterator<Q1, 1, device::device_pointer<Q2>> last,
                                                   multi::array_iterator<T, 1, T*> dest)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  assert(stride(first) == stride(last));
  if (std::distance(first, last) == 0)
    return dest;
  using qmcplusplus::afqmc::to_address;
  arch::memcopy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                  sizeof(T), std::distance(first, last), arch::memcopyD2H);
  return dest + std::distance(first, last);
}

template<class Alloc, class T, class Q1, class Q2>
multi::array_iterator<T, 1, T*> alloc_uninitialized_copy(Alloc& a,
                                                         multi::array_iterator<Q1, 1, device::device_pointer<Q2>> first,
                                                         multi::array_iterator<Q1, 1, device::device_pointer<Q2>> last,
                                                         multi::array_iterator<T, 1, T*> dest)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  assert(stride(first) == stride(last));
  if (std::distance(first, last) == 0)
    return dest;
  using qmcplusplus::afqmc::to_address;
  arch::memcopy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                  sizeof(T), std::distance(first, last), arch::memcopyD2H);
  return dest + std::distance(first, last);
}

/*
template<class Alloc, class T, class Q, typename Size>
multi::array_iterator<T, 1, device::device_pointer<T>> uninitialized_copy_n( 
                           Alloc &a,
                           multi::array_iterator<Q, 1, device::device_pointer<Q>> first,
                           Size N,
                           multi::array_iterator<T, 1, device::device_pointer<T>> dest ){
  static_assert(std::is_same<typename std::decay<Q>::type,T>::value,"Wrong dispatch.\n");
  if(N==0) return dest;
  if(cudaSuccess != cudaMemcpy2D(to_address(base(dest)),sizeof(T)*stride(dest),
                                 to_address(base(first)),sizeof(T)*stride(first),
                                 sizeof(T),N,cudaMemcpyDeviceToDevice))
      throw std::runtime_error("Error: cudaMemcpy2D returned error code.");
  return dest+N;
}
*/

template<class Alloc, class T, class Q1, class Q2, typename Size>
multi::array_iterator<T, 1, T*> uninitialized_copy_n(Alloc& a,
                                                     multi::array_iterator<Q1, 1, device::device_pointer<Q2>> first,
                                                     Size n,
                                                     multi::array_iterator<T, 1, T*> dest)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  if (n == 0)
    return dest;
  using qmcplusplus::afqmc::to_address;
  arch::memcopy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                  sizeof(T), n);
  return dest + n;
}

template<class Alloc, class T, class ForwardIt, typename Size>
multi::array_iterator<T, 1, device::device_pointer<T>> uninitialized_copy_n(
    Alloc& a,
    ForwardIt first,
    Size n,
    multi::array_iterator<T, 1, device::device_pointer<T>> dest)
{
  if (n == 0)
    return dest;
  using qmcplusplus::afqmc::to_address;
  arch::memcopy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                  sizeof(T), n);
  return dest + n;
}

template<class Alloc, typename T, typename Size>
multi::array_iterator<T, 1, device::device_pointer<T>> uninitialized_default_construct_n(
    Alloc& a,
    multi::array_iterator<T, 1, device::device_pointer<T>> first,
    Size n)
{
  return uninitialized_fill_n(first, n, T());
}

template<class Alloc, typename T>
multi::array_iterator<T, 1, device::device_pointer<T>> uninitialized_default_construct(
    Alloc& a,
    multi::array_iterator<T, 1, device::device_pointer<T>> first,
    multi::array_iterator<T, 1, device::device_pointer<T>> last)
{
  return uninitialized_fill_n(a, first, std::distance(first, last), T());
}

template<class Alloc, typename T, typename Size>
multi::array_iterator<T, 1, device::device_pointer<T>> uninitialized_value_construct_n(
    Alloc& a,
    multi::array_iterator<T, 1, device::device_pointer<T>> first,
    Size n)
{
  return uninitialized_fill_n(first, n, T());
}

template<class Alloc, typename T>
multi::array_iterator<T, 1, device::device_pointer<T>> uninitialized_value_construct(
    Alloc& a,
    multi::array_iterator<T, 1, device::device_pointer<T>> first,
    multi::array_iterator<T, 1, device::device_pointer<T>> last)
{
  return uninitialized_fill_n(a, first, std::distance(first, last), T());
}

template<class Alloc, class T, class Q1, class Q2, typename Size>
multi::array_iterator<T, 1, T*> alloc_uninitialized_copy_n(
    Alloc& a,
    multi::array_iterator<Q1, 1, device::device_pointer<Q2>> first,
    Size n,
    multi::array_iterator<T, 1, T*> dest)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  if (n == 0)
    return dest;
  using qmcplusplus::afqmc::to_address;
  arch::memcopy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                  sizeof(T), n);
  return dest + n;
}

template<class Alloc, class T, class ForwardIt, typename Size>
multi::array_iterator<T, 1, device::device_pointer<T>> alloc_uninitialized_copy_n(
    Alloc& a,
    ForwardIt first,
    Size n,
    multi::array_iterator<T, 1, device::device_pointer<T>> dest)
{
  if (n == 0)
    return dest;
  using qmcplusplus::afqmc::to_address;
  arch::memcopy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                  sizeof(T), n);
  return dest + n;
}

template<class Alloc, class T, class Q1, class Q2, typename Size>
multi::array_iterator<T, 1, T*> alloc_uninitialized_move_n(
    Alloc& a,
    multi::array_iterator<Q1, 1, device::device_pointer<Q2>> first,
    Size n,
    multi::array_iterator<T, 1, T*> dest)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  if (n == 0)
    return dest;
  using qmcplusplus::afqmc::to_address;
  arch::memcopy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                  sizeof(T), n);
  return dest + n;
}

template<class Alloc, class T, class ForwardIt, typename Size>
multi::array_iterator<T, 1, device::device_pointer<T>> alloc_uninitialized_move_n(
    Alloc& a,
    ForwardIt first,
    Size n,
    multi::array_iterator<T, 1, device::device_pointer<T>> dest)
{
  if (n == 0)
    return dest;
  using qmcplusplus::afqmc::to_address;
  arch::memcopy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                  sizeof(T), n);
  return dest + n;
}

template<class Alloc, typename T, typename Size>
multi::array_iterator<T, 1, device::device_pointer<T>> alloc_uninitialized_default_construct_n(
    Alloc& a,
    multi::array_iterator<T, 1, device::device_pointer<T>> first,
    Size n)
{
  return uninitialized_fill_n(first, n, T());
}

template<class Alloc, typename T>
multi::array_iterator<T, 1, device::device_pointer<T>> alloc_uninitialized_default_construct(
    Alloc& a,
    multi::array_iterator<T, 1, device::device_pointer<T>> first,
    multi::array_iterator<T, 1, device::device_pointer<T>> last)
{
  return uninitialized_fill_n(a, first, std::distance(first, last), T());
}

template<class Alloc, typename T, typename Size>
multi::array_iterator<T, 1, device::device_pointer<T>> alloc_uninitialized_value_construct_n(
    Alloc& a,
    multi::array_iterator<T, 1, device::device_pointer<T>> first,
    Size n)
{
  return uninitialized_fill_n(first, n, T());
}

template<class Alloc, typename T>
multi::array_iterator<T, 1, device::device_pointer<T>> alloc_uninitialized_value_construct(
    Alloc& a,
    multi::array_iterator<T, 1, device::device_pointer<T>> first,
    multi::array_iterator<T, 1, device::device_pointer<T>> last)
{
  return uninitialized_fill_n(a, first, std::distance(first, last), T());
}


} // namespace multi
} // namespace boost


#endif
