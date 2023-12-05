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

#ifndef AFQMC_CUDA_UM_POINTERS_HPP
#define AFQMC_CUDA_UM_POINTERS_HPP

#include <cassert>
#include "Numerics/detail/cuda_utilities.hpp"

namespace cuda
{
template<class T>
struct cuda_um_ptr
{
  using value_type             = T;
  using const_value_type       = T const;
  using pointer                = T*;
  using const_pointer          = T const*;
  static const int memory_type = MANAGED_MEMORY_POINTER_TYPE;
  T* impl_;
  cuda::gpu_handles handles;
  cuda_um_ptr() = default;
  cuda_um_ptr(T* impl__, cuda::gpu_handles handles_ = cuda::gpu_handles{}) : impl_(impl__), handles(handles_) {}
  // eventually check if memory types and blas types are convertible, e.g. CPU_MEMORY to CPU_OUTOFCARD
  template<typename Q, typename = typename std::enable_if_t<cuda_um_ptr<Q>::memory_type == memory_type>>
  cuda_um_ptr(cuda_um_ptr<Q> const& ptr) : impl_(ptr.impl_), handles(ptr.handles)
  {}
  T& operator*() const { return *impl_; }
  T& operator[](std::ptrdiff_t n) const { return impl_[n]; }
  T* operator->() const { return impl_; }
  explicit operator bool() const { return (impl_ != nullptr); }
  auto operator+(std::ptrdiff_t n) const { return cuda_um_ptr{impl_ + n, handles}; }
  std::ptrdiff_t operator-(cuda_um_ptr other) const { return std::ptrdiff_t(impl_ - other.impl_); }
  cuda_um_ptr& operator++()
  {
    ++impl_;
    return *this;
  }
  cuda_um_ptr& operator--()
  {
    --impl_;
    return *this;
  }
  cuda_um_ptr& operator+=(std::ptrdiff_t d)
  {
    impl_ += d;
    return *this;
  }
  cuda_um_ptr& operator-=(std::ptrdiff_t d)
  {
    impl_ -= d;
    return *this;
  }
  bool operator==(cuda_um_ptr const& other) const { return impl_ == other.impl_; }
  bool operator!=(cuda_um_ptr const& other) const { return not(*this == other); }
  bool operator<=(cuda_um_ptr<T> const& other) const { return impl_ <= other.impl_; }
  T* to_address() const { return impl_; }
  friend decltype(auto) to_address(cuda_um_ptr const& self) { return self.to_address(); }
  template<class Q>
  friend cuda_um_ptr<Q> pointer_cast(cuda_um_ptr<T> const& self)
  {
    return cuda_um_ptr<Q>{reinterpret_cast<Q*>(self.impl_), self.handles};
  }
};


template<class T>
struct cuda_um_allocator
{
  template<class U>
  struct rebind
  {
    using other = cuda_um_allocator<U>;
  };
  using value_type       = T;
  using const_value_type = T const;
  using pointer          = cuda_um_ptr<T>;
  using const_pointer    = cuda_um_ptr<T const>;
  using reference        = T&;
  using const_reference  = T const&;
  using size_type        = std::size_t;
  using difference_type  = std::ptrdiff_t;

  cuda::gpu_handles handles_;
  cuda_um_allocator(cuda::gpu_handles handles__) : handles_(handles__) {}
  cuda_um_allocator()  = delete;
  ~cuda_um_allocator() = default;
  cuda_um_allocator(cuda_um_allocator const& other) : handles_(other.handles_) {}
  template<class U>
  cuda_um_allocator(cuda_um_allocator<U> const& other) : handles_(other.handles_)
  {}

  cuda_um_ptr<T> allocate(size_type n, const void* hint = 0)
  {
    if (n == 0)
      return cuda_um_ptr<T>{};
    T* p;
    cudaMallocManaged((void**)&p, n * sizeof(T), cudaMemAttachGlobal);
    return cuda_um_ptr<T>{p, handles_};
  }
  void deallocate(cuda_um_ptr<T> ptr, size_type) { cudaFree(ptr.impl_); }
  bool operator==(cuda_um_allocator const& other) const { return (handles_ == other.handles_); }
  bool operator!=(cuda_um_allocator const& other) const { return not(other == *this); }
  template<class U, class... Args>
  void construct(U* p, Args&&... args)
  {
    ::new ((void*)p) U(std::forward<Args>(args)...);
  }
  template<class U>
  void destroy(U* p)
  {
    p->~U();
  }
};

template<class T, class F>
F for_each(cuda_um_ptr<T> first, cuda_um_ptr<T> last, F f)
{
  if (first == last)
    return f;
  return std::for_each(to_address(first), to_address(last), f);
}

template<typename T, typename Size, typename... Args>
cuda_um_ptr<T> fill_n(cuda_um_ptr<T> first, Size n, Args&&... args)
{
  if (n == 0)
    return first;
  std::fill_n(to_address(first), n, std::forward<Args>(args)...);
  return first + n;
}

template<typename T, typename Size, typename... Args>
cuda_um_ptr<T> uninitialized_fill_n(cuda_um_ptr<T> first, Size n, Args&&... args)
{
  if (n == 0)
    return first;
  std::uninitialized_fill_n(to_address(first), n, std::forward<Args>(args)...);
  return first + n;
}

template<typename T, typename Size>
cuda_um_ptr<T> destroy_n(cuda_um_ptr<T> first, Size n)
{
  auto first_ptr = to_address(first);
  for (; n > 0; (void)++first_ptr, --n)
    first->~T();
  return first + n;
}

} // namespace cuda

#endif
