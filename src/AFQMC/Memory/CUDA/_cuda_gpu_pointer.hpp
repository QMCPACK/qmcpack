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

#ifndef AFQMC_CUDA_GPU_POINTERS_HPP
#define AFQMC_CUDA_GPU_POINTERS_HPP

#include <functional>
#include "Configuration.h"
#include <cassert>
#include <cuda_runtime.h>
//#include "cublas_v2.h"
//#include "cublasXt.h"
#include "AFQMC/Memory/CUDA/cuda_utilities.h"
#include "AFQMC/Numerics/detail/CUDA/Kernels/fill_n.cuh"
//#include "AFQMC/Numerics/detail/CUDA/Kernels/uninitialized_fill_n.cuh"
//#include "AFQMC/Numerics/detail/CUDA/Kernels/uninitialized_copy_n.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/copy_n_cast.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/print.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/reference_operations.cuh"

#include "multi/array_ref.hpp"
#include "mpi3/shared_window.hpp"
#include "AFQMC/Memory/raw_pointers.hpp"

namespace qmc_cuda
{
using qmcplusplus::afqmc::to_address;

template<class T>
struct cuda_gpu_allocator;
template<class T>
struct device_pointer;

// no const for now
template<class T>
struct cuda_gpu_reference
{
public:
  using value_type       = T;
  using decay_value_type = typename std::decay<T>::type;
  using pointer          = device_pointer<T>;

  // must construct through a gpu_ptr for now, to keep some sort of control/safety
  cuda_gpu_reference(pointer const& gpu_ptr) : impl_(to_address(gpu_ptr)) {}
  cuda_gpu_reference(cuda_gpu_reference<T> const& gpu_ref) = default;

  // assignment
  cuda_gpu_reference& operator=(cuda_gpu_reference const& x)
  {
    if (cudaSuccess != cudaMemcpy(impl_, x.impl_, sizeof(T), cudaMemcpyDefault))
      throw std::runtime_error("Error: cudaMemcpy returned error code.");
    return *this;
  }

  cuda_gpu_reference& operator=(value_type const& x)
  {
    if (cudaSuccess != cudaMemcpy(impl_, std::addressof(x), sizeof(T), cudaMemcpyDefault))
      throw std::runtime_error("Error: cudaMemcpy returned error code.");
    return *this;
  }

  operator value_type() const { return this->val(); }

  // try getting rid of this
  operator value_type&()
  {
    if (cudaSuccess != cudaMemcpy(std::addressof(host_impl_), impl_, sizeof(T), cudaMemcpyDefault))
      throw std::runtime_error("Error: cudaMemcpy returned error code.");
    return host_impl_;
  }

  pointer operator&() const { return pointer{impl_}; }

  value_type* addressof() { return impl_; }
  value_type const* addressof() const { return impl_; }
  friend auto addressof(cuda_gpu_reference const& rhs) { return rhs.addressof(); }
  friend auto addressof(cuda_gpu_reference&& rhs) { return rhs.addressof(); }
  friend auto addressof(cuda_gpu_reference& rhs) { return rhs.addressof(); }

  void swap(cuda_gpu_reference& other) { std::swap(impl_, other.impl_); }

  cuda_gpu_reference& operator++(void)
  {
    kernels::op_plus(impl_, T(1));
    return *this;
  }
  cuda_gpu_reference& operator--(void)
  {
    kernels::op_minus(impl_, T(1));
    return *this;
  }

  value_type operator++(int)
  {
    value_type res;
    if (cudaSuccess != cudaMemcpy(std::addressof(res), impl_, sizeof(T), cudaMemcpyDefault))
      throw std::runtime_error("Error: cudaMemcpy returned error code.");
    return res++;
  }
  value_type operator--(int)
  {
    value_type res;
    if (cudaSuccess != cudaMemcpy(std::addressof(res), impl_, sizeof(T), cudaMemcpyDefault))
      throw std::runtime_error("Error: cudaMemcpy returned error code.");
    return res--;
  }

  cuda_gpu_reference& operator+=(const value_type& rhs)
  {
    kernels::op_plus(impl_, rhs);
    return *this;
  }
  cuda_gpu_reference& operator-=(const value_type& rhs)
  {
    kernels::op_minus(impl_, rhs);
    return *this;
  }
  cuda_gpu_reference& operator*=(const value_type& rhs)
  {
    kernels::op_times(impl_, rhs);
    return *this;
  }
  cuda_gpu_reference& operator/=(const value_type& rhs)
  {
    kernels::op_div(impl_, rhs);
    return *this;
  }

  cuda_gpu_reference& operator+=(cuda_gpu_reference const& rhs)
  {
    kernels::op_plus(impl_, rhs.val());
    return *this;
  }
  cuda_gpu_reference& operator-=(cuda_gpu_reference const& rhs)
  {
    kernels::op_minus(impl_, rhs.val());
    return *this;
  }
  cuda_gpu_reference& operator*=(cuda_gpu_reference const& rhs)
  {
    kernels::op_times(impl_, rhs.val());
    return *this;
  }
  cuda_gpu_reference& operator/=(cuda_gpu_reference const& rhs)
  {
    kernels::op_div(impl_, rhs.val());
    return *this;
  }

  friend value_type& operator+=(value_type& lhs, cuda_gpu_reference const& rhs)
  {
    lhs += rhs.val();
    return lhs;
  }
  friend value_type& operator-=(value_type& lhs, cuda_gpu_reference const& rhs)
  {
    lhs -= rhs.val();
    return lhs;
  }
  friend value_type& operator*=(value_type& lhs, cuda_gpu_reference const& rhs)
  {
    lhs *= rhs.val();
    return lhs;
  }
  friend value_type& operator/=(value_type& lhs, cuda_gpu_reference const& rhs)
  {
    lhs /= rhs.val();
    return lhs;
  }

  value_type operator+(value_type const& rhs) const { return this->val() + rhs; }
  value_type operator-(value_type const& rhs) const { return this->val() - rhs; }
  value_type operator/(value_type const& rhs) const { return this->val() / rhs; }
  value_type operator*(value_type const& rhs) const { return this->val() * rhs; }

  value_type operator+(cuda_gpu_reference const& rhs) const { return this->val() + rhs.val(); }
  value_type operator-(cuda_gpu_reference const& rhs) const { return this->val() - rhs.val(); }
  value_type operator*(cuda_gpu_reference const& rhs) const { return this->val() * rhs.val(); }
  value_type operator/(cuda_gpu_reference const& rhs) const { return this->val() / rhs.val(); }

  friend value_type operator+(value_type lhs, cuda_gpu_reference const& rhs) { return lhs + rhs.val(); }
  friend value_type operator-(value_type lhs, cuda_gpu_reference const& rhs) { return lhs - rhs.val(); }
  friend value_type operator*(value_type lhs, cuda_gpu_reference const& rhs) { return lhs * rhs.val(); }
  friend value_type operator/(value_type lhs, cuda_gpu_reference const& rhs) { return lhs / rhs.val(); }

  bool operator==(value_type const& rhs) const { return this->val() == rhs; }
  bool operator!=(value_type const& rhs) const { return this->val() != rhs; }
  bool operator>(value_type const& rhs) const { return this->val() > rhs; }
  bool operator<(value_type const& rhs) const { return this->val() < rhs; }
  bool operator>=(value_type const& rhs) const { return this->val() >= rhs; }
  bool operator<=(value_type const& rhs) const { return this->val() <= rhs; }

  bool operator==(cuda_gpu_reference const& rhs) const { return this->val() == rhs.val(); }
  bool operator!=(cuda_gpu_reference const& rhs) const { return this->val() != rhs.val(); }
  bool operator>(cuda_gpu_reference const& rhs) const { return this->val() > rhs.val(); }
  bool operator<(cuda_gpu_reference const& rhs) const { return this->val() < rhs.val(); }
  bool operator>=(cuda_gpu_reference const& rhs) const { return this->val() >= rhs.val(); }
  bool operator<=(cuda_gpu_reference const& rhs) const { return this->val() <= rhs.val(); }

  friend bool operator==(value_type const& lhs, cuda_gpu_reference const& rhs) { return lhs == rhs.val(); }
  friend bool operator!=(value_type const& lhs, cuda_gpu_reference const& rhs) { return lhs != rhs.val(); }
  friend bool operator>(value_type const& lhs, cuda_gpu_reference const& rhs) { return lhs > rhs.val(); }
  friend bool operator<(value_type const& lhs, cuda_gpu_reference const& rhs) { return lhs < rhs.val(); }
  friend bool operator>=(value_type const& lhs, cuda_gpu_reference const& rhs) { return lhs >= rhs.val(); }
  friend bool operator<=(value_type const& lhs, cuda_gpu_reference const& rhs) { return lhs <= rhs.val(); }

  friend std::ostream& operator<<(std::ostream& os, cuda_gpu_reference const& obj)
  {
    os << obj.val();
    return os;
  }

  friend std::istream& operator<<(std::istream& is, cuda_gpu_reference& obj)
  {
    value_type val;
    is >> val;
    obj = val;
    return is;
  }

  /*
  cuda_gpu_reference&  operator<<= (const T &rhs)
 
  cuda_gpu_reference&  operator>>= (const T &rhs)
 
  cuda_gpu_reference&  operator&= (const T &rhs)
 
  cuda_gpu_reference&  operator|= (const T &rhs)
 
  cuda_gpu_reference&  operator^= (const T &rhs)
*/

private:
  value_type* impl_;
  //value_type host_impl_;
  decay_value_type host_impl_;

  value_type val() const
  {
    decay_value_type res;
    if (cudaSuccess != cudaMemcpy(std::addressof(res), impl_, sizeof(value_type), cudaMemcpyDefault))
      throw std::runtime_error("Error: cudaMemcpy returned error code.");
    return value_type(res);
  }
};


struct base_device_pointer
{
  static gpu_handles handles;
};

// this class is not safe, since it allows construction of a gpu_ptr from a raw ptr
// which might not be in gpu memory. Fix this!!!
template<class T>
struct device_pointer : base_device_pointer
{
  using difference_type  = std::ptrdiff_t;
  using value_type       = T;
  using element_type     = T;
  using const_value_type = T const;
  using pointer          = T*;
  using const_pointer    = T const*;
  // this is wrong!!! but no synthetic references yet!!!
  using reference              = cuda_gpu_reference<T>;
  using const_reference        = cuda_gpu_reference<T> const;
  using iterator_category      = std::random_access_iterator_tag;
  static const int memory_type = GPU_MEMORY_POINTER_TYPE;
  using default_allocator_type = cuda_gpu_allocator<T>;
  friend class cuda_gpu_allocator<T>;
  friend class device_pointer<typename std::decay<T>::type>;
  default_allocator_type default_allocator() const { return cuda_gpu_allocator<T>{}; };
  device_pointer() = default;
  device_pointer(std::nullptr_t) : impl_(nullptr) {}
  // eventually check if memory types and blas types are convertible, e.g. CPU_MEMORY to CPU_OUTOFCARD
  template<typename Q>
  device_pointer(device_pointer<Q> const& ptr) : impl_(ptr.impl_)
  {}
  reference operator*() const { return reference(device_pointer{impl_}); }
  reference operator[](std::ptrdiff_t n) const { return reference(device_pointer{impl_ + n}); }
  T* operator->() const { return impl_; }
  explicit operator bool() const { return (impl_ != nullptr); }
  //  operator device_pointer<T const>() const{return device_pointer<T const>{impl_}; }
  auto operator+(std::ptrdiff_t n) const { return device_pointer{impl_ + n}; }
  std::ptrdiff_t operator-(device_pointer other) const { return std::ptrdiff_t(impl_ - other.impl_); }
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

protected:
  device_pointer(T* impl__) : impl_(impl__) {}
};

//static size_t TotalGPUAlloc=0;

/*
 * Incomplete!!! Need to fix construct and destroy
 */
template<class T>
struct cuda_gpu_allocator
{
  template<class U>
  struct rebind
  {
    using other = cuda_gpu_allocator<U>;
  };
  using element_type     = T;
  using value_type       = T;
  using const_value_type = T const;
  using pointer          = device_pointer<T>;
  using const_pointer    = device_pointer<T const>;
  using reference        = T&;
  using const_reference  = T const&;
  using size_type        = std::size_t;
  using difference_type  = std::ptrdiff_t;

  cuda_gpu_allocator()                                = default;
  ~cuda_gpu_allocator()                               = default;
  cuda_gpu_allocator(cuda_gpu_allocator const& other) = default;
  template<class U>
  cuda_gpu_allocator(cuda_gpu_allocator<U> const& other)
  {}

  device_pointer<T> allocate(size_type n, const void* hint = 0)
  {
    if (n == 0)
      return device_pointer<T>{};
    T* p;
    if (cudaSuccess != cudaMalloc((void**)&p, n * sizeof(T)))
    {
      std::cerr << " Error allocating " << n * sizeof(T) / 1024.0 / 1024.0 << " MBs on GPU." << std::endl;
      throw std::runtime_error("Error: cudaMalloc returned error code.");
    }
    device_pointer<T> res;
    res.impl_ = p;
    return res;
    //    return device_pointer<T>{p};
  }
  void deallocate(device_pointer<T> ptr, size_type) { cudaFree(ptr.impl_); }
  bool operator==(cuda_gpu_allocator const& other) const { return true; }
  bool operator!=(cuda_gpu_allocator const& other) const { return false; }
  template<class U, class... Args>
  void construct(U* p, Args&&... args)
  {}
  //{
  //    ::new((void*)p) U(std::forward<Args>(args)...);
  //  }
  template<class U>
  void destroy(U* p)
  {}
  //  {
  //    p->~U();
  //  }
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
template<typename T, typename Size>
device_pointer<T> copy_n(device_pointer<T> const A, Size n, device_pointer<T> B)
{
  if (cudaSuccess != cudaMemcpy(to_address(B), to_address(A), n * sizeof(T), cudaMemcpyDefault))
    throw std::runtime_error("Error: cudaMemcpy returned error code.");
  return B + n;
}

template<typename T, typename Size>
device_pointer<T> copy_n(T* const A, Size n, device_pointer<T> B)
{
  if (cudaSuccess != cudaMemcpy(to_address(B), A, n * sizeof(T), cudaMemcpyDefault))
    throw std::runtime_error("Error: cudaMemcpy returned error code.");
  return B + n;
}

template<typename T, typename Size>
T* copy_n(device_pointer<T> const A, Size n, T* B)
{
  if (cudaSuccess != cudaMemcpy(B, to_address(A), n * sizeof(T), cudaMemcpyDefault))
    throw std::runtime_error("Error: cudaMemcpy returned error code.");
  return B + n;
}

/**************** copy *****************/
template<typename T>
device_pointer<T> copy(device_pointer<T> const Abeg, device_pointer<T> const Aend, device_pointer<T> B)
{
  return copy_n(Abeg, std::distance(Abeg, Aend), B);
}

template<typename T>
device_pointer<T> copy(T* const Abeg, T* const Aend, device_pointer<T> B)
{
  return copy_n(Abeg, std::distance(Abeg, Aend), B);
}

template<typename T>
T* copy(device_pointer<T> const Abeg, device_pointer<T> const Aend, T* B)
{
  return copy_n(Abeg, std::distance(Abeg, Aend), B);
}

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
  throw std::runtime_error(" Error: copy_n_cast(gpu_ptr,n,T*) is disabled.");
  return B + n;
}

template<typename T, typename Q, typename Size>
Q* copy_n_cast(device_pointer<T> const A, Size n, Q* B)
{
  throw std::runtime_error(" Error: copy_n_cast(gpu_ptr,n,T*) is disabled.");
  return B + n;
}

/**************** fill_n *****************/
//template<typename T, typename Size, typename... Args>
//device_pointer<T> fill_n(device_pointer<T> first, Size n, Args&&...args){
//  if(n == 0) return first;
//  kernels::fill_n(to_address(first), n, std::forward<Args>(args)...);
//  return first + n;
//}

template<typename T, typename Size>
device_pointer<T> fill_n(device_pointer<T> first, Size n, T const& val)
{
  if (n == 0)
    return first;
  kernels::fill_n(to_address(first), n, val);
  return first + n;
}

template<typename T>
device_pointer<T> fill(device_pointer<T> first, device_pointer<T> last, T const& val)
{
  return fill_n(first, std::distance(first, last), val);
}


/**************** uninitialized_fill_n *****************/

template<typename T, typename Size>
device_pointer<T> uninitialized_fill_n(device_pointer<T> first, Size n, T const& val)
{
  if (n == 0)
    return first;
  //kernels::uninitialized_fill_n(to_address(first), n, val);
  kernels::fill_n(to_address(first), n, val);
  return first + n;
}

template<typename T>
device_pointer<T> uninitialized_fill(device_pointer<T> first, device_pointer<T> last, T const& val)
{
  return uninitialized_fill_n(first, std::distance(first, last), val);
}

template<class Alloc, typename T, typename Size>
device_pointer<T> uninitialized_fill_n(Alloc& a, device_pointer<T> first, Size n, T const& val)
{
  if (n == 0)
    return first;
  //kernels::uninitialized_fill_n(to_address(first), n, val);
  kernels::fill_n(to_address(first), n, val);
  return first + n;
}

template<class Alloc, typename T>
device_pointer<T> uninitialized_fill(Alloc& a, device_pointer<T> first, device_pointer<T> last, T const& val)
{
  return uninitialized_fill_n(a, first, std::distance(first, last), val);
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

template<typename T, typename Size>
device_pointer<T> uninitialized_copy_n(T* A, Size n, device_pointer<T> B)
{
  return copy_n(A, n, B);
}

template<typename T, typename Size>
T* uninitialized_copy_n(device_pointer<T> A, Size n, T* B)
{
  return copy_n(A, n, B);
}

template<class Alloc, typename T, typename Size>
device_pointer<T> uninitialized_copy_n(Alloc& a, device_pointer<T> A, Size n, device_pointer<T> B)
{
  return copy_n(A, n, B);
}

template<class Alloc, typename T, typename Size>
device_pointer<T> uninitialized_copy_n(Alloc& a, T* A, Size n, device_pointer<T> B)
{
  return copy_n(A, n, B);
}

template<class Alloc, typename T, typename Size>
T* uninitialized_copy_n(Alloc& a, device_pointer<T> A, Size n, T* B)
{
  return copy_n(A, n, B);
}

/**************** uninitialized_copy *****************/
template<typename T>
device_pointer<T> uninitialized_copy(device_pointer<T> const Abeg, device_pointer<T> const Aend, device_pointer<T> B)
{
  return copy_n(Abeg, std::distance(Abeg, Aend), B);
}

template<typename T>
device_pointer<T> uninitialized_copy(T* const Abeg, T* const Aend, device_pointer<T> B)
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

/**************** print *****************/
template<typename T>
void print(std::string str, device_pointer<T> p, int n)
{
  kernels::print(str, to_address(p), n);
}

template<typename T>
void fill2D(int n, int m, qmc_cuda::device_pointer<T> first, int lda, T const& val)
{
  assert(lda >= m);
  kernels::fill2D_n(n, m, to_address(first), lda, val);
}

} // namespace qmc_cuda

namespace boost
{
namespace multi
{
//using qmcplusplus::afqmc::to_address;

/**************** copy *****************/
// Can always call cudaMemcopy2D like you do in the blas backend

template<class Alloc, typename T, typename Size>
multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> uninitialized_fill_n(
    Alloc& a,
    multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> first,
    Size n,
    T const& val)
{
  if (n == 0)
    return first;
  kernels::fill_n(to_address(base(first)), n, stride(first), val);
  return first + n;
}

template<class Alloc, typename T, typename Size>
multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> uninitialized_fill(
    Alloc& a,
    multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> first,
    multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> last,
    T const& val)
{
  assert(stride(first) == stride(last));
  if (std::distance(first, last) == 0)
    return first;
  kernels::fill_n(to_address(base(first)), std::distance(first, last), stride(first), val);
  return first + std::distance(first, last);
}

template<class T, class Q1, class Q2>
multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> copy(
    multi::array_iterator<Q1, 1, qmc_cuda::device_pointer<Q2>> first,
    multi::array_iterator<Q1, 1, qmc_cuda::device_pointer<Q2>> last,
    multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> dest)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  assert(stride(first) == stride(last));
  if (std::distance(first, last) == 0)
    return dest;
  if (cudaSuccess !=
      cudaMemcpy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                   sizeof(T), std::distance(first, last), cudaMemcpyDeviceToDevice))
    throw std::runtime_error("Error: cudaMemcpy2D returned error code.");
  return dest + std::distance(first, last);
}

template<class T, class ForwardIt>
multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> copy(
    ForwardIt first,
    ForwardIt last,
    multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> dest)
{
  assert(stride(first) == stride(last));
  if (std::distance(first, last) == 0)
    return dest;
  if (cudaSuccess !=
      cudaMemcpy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                   sizeof(T), std::distance(first, last), cudaMemcpyDefault))
    throw std::runtime_error("Error: cudaMemcpy2D returned error code.");
  return dest + std::distance(first, last);
}

template<class T, class Q1, class Q2>
multi::array_iterator<T, 1, T*> copy(multi::array_iterator<Q1, 1, qmc_cuda::device_pointer<Q2>> first,
                                     multi::array_iterator<Q1, 1, qmc_cuda::device_pointer<Q2>> last,
                                     multi::array_iterator<T, 1, T*> dest)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  assert(stride(first) == stride(last));
  if (std::distance(first, last) == 0)
    return dest;
  if (cudaSuccess !=
      cudaMemcpy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                   sizeof(T), std::distance(first, last), cudaMemcpyDeviceToHost))
    throw std::runtime_error("Error: cudaMemcpy2D returned error code.");
  return dest + std::distance(first, last);
}

template<class T, class Q1, class Q2, typename Size>
multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> copy_n(
    multi::array_iterator<Q1, 1, qmc_cuda::device_pointer<Q2>> first,
    Size N,
    multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> dest)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  if (N == 0)
    return dest;
  if (cudaSuccess !=
      cudaMemcpy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                   sizeof(T), N, cudaMemcpyDeviceToDevice))
    throw std::runtime_error("Error: cudaMemcpy2D returned error code.");
  return dest + N;
}

template<class T, class ForwardIt, typename Size>
multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> copy_n(
    ForwardIt first,
    Size n,
    multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> dest)
{
  if (n == 0)
    return dest;
  if (cudaSuccess !=
      cudaMemcpy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                   sizeof(T), n, cudaMemcpyDefault))
    throw std::runtime_error("Error: cudaMemcpy2D returned error code.");
  return dest + n;
}

template<class T, class Q1, class Q2, typename Size>
multi::array_iterator<T, 1, T*> copy_n(multi::array_iterator<Q1, 1, qmc_cuda::device_pointer<Q2>> first,
                                       Size N,
                                       multi::array_iterator<T, 1, T*> dest)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  if (N == 0)
    return dest;
  if (cudaSuccess !=
      cudaMemcpy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                   sizeof(T), N, cudaMemcpyDeviceToHost))
    throw std::runtime_error("Error: cudaMemcpy2D returned error code.");
  return dest + N;
}

template<class Alloc, class T, class ForwardIt>
multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> uninitialized_copy(
    Alloc& a,
    ForwardIt first,
    ForwardIt last,
    multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> dest)
{
  assert(stride(first) == stride(last));
  if (std::distance(first, last) == 0)
    return dest;
  if (cudaSuccess !=
      cudaMemcpy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                   sizeof(T), std::distance(first, last), cudaMemcpyDefault))
    throw std::runtime_error("Error: cudaMemcpy2D returned error code.");
  return dest + std::distance(first, last);
}

/*
template<class Alloc, class T, class Q>
multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> uninitialized_copy( 
                         Alloc &a,
                         multi::array_iterator<Q, 1, qmc_cuda::device_pointer<Q>> first,
                         multi::array_iterator<Q, 1, qmc_cuda::device_pointer<Q>> last,
                         multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> dest ){
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
multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> uninitialized_copy(
                         Alloc &a,
                         multi::array_iterator<Q, 1, boost::mpi3::intranode::array_ptr<Q>> first,
                         multi::array_iterator<Q, 1, boost::mpi3::intranode::array_ptr<Q>> last,
                         multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> dest ){
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
multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> uninitialized_copy(
                         Alloc &a,
                         multi::array_iterator<Q1, 1, Q2*> first,
                         multi::array_iterator<Q1, 1, Q2*> last,
                         multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> dest ){
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
                                                   multi::array_iterator<Q1, 1, qmc_cuda::device_pointer<Q2>> first,
                                                   multi::array_iterator<Q1, 1, qmc_cuda::device_pointer<Q2>> last,
                                                   multi::array_iterator<T, 1, T*> dest)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  assert(stride(first) == stride(last));
  if (std::distance(first, last) == 0)
    return dest;
  using qmcplusplus::afqmc::to_address;
  if (cudaSuccess !=
      cudaMemcpy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                   sizeof(T), std::distance(first, last), cudaMemcpyDeviceToHost))
    throw std::runtime_error("Error: cudaMemcpy2D returned error code.");
  return dest + std::distance(first, last);
}

/*
template<class Alloc, class T, class Q, typename Size>
multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> uninitialized_copy_n( 
                           Alloc &a,
                           multi::array_iterator<Q, 1, qmc_cuda::device_pointer<Q>> first,
                           Size N,
                           multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> dest ){
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
                                                     multi::array_iterator<Q1, 1, qmc_cuda::device_pointer<Q2>> first,
                                                     Size n,
                                                     multi::array_iterator<T, 1, T*> dest)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  if (n == 0)
    return dest;
  if (cudaSuccess !=
      cudaMemcpy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                   sizeof(T), n, cudaMemcpyDefault))
    throw std::runtime_error("Error: cudaMemcpy2D returned error code.");
  return dest + n;
}

template<class Alloc, class T, class ForwardIt, typename Size>
multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> uninitialized_copy_n(
    Alloc& a,
    ForwardIt first,
    Size n,
    multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> dest)
{
  if (n == 0)
    return dest;
  if (cudaSuccess !=
      cudaMemcpy2D(to_address(base(dest)), sizeof(T) * stride(dest), to_address(base(first)), sizeof(T) * stride(first),
                   sizeof(T), n, cudaMemcpyDefault))
    throw std::runtime_error("Error: cudaMemcpy2D returned error code.");
  return dest + n;
}

template<class Alloc, typename T, typename Size>
multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> uninitialized_default_construct_n(
    Alloc& a,
    multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> first,
    Size n)
{
  return uninitialized_fill_n(first, n, T());
}

template<class Alloc, typename T>
multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> uninitialized_default_construct(
    Alloc& a,
    multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> first,
    multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> last)
{
  return uninitialized_fill_n(a, first, std::distance(first, last), T());
}

template<class Alloc, typename T, typename Size>
multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> uninitialized_value_construct_n(
    Alloc& a,
    multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> first,
    Size n)
{
  return uninitialized_fill_n(first, n, T());
}

template<class Alloc, typename T>
multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> uninitialized_value_construct(
    Alloc& a,
    multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> first,
    multi::array_iterator<T, 1, qmc_cuda::device_pointer<T>> last)
{
  return uninitialized_fill_n(a, first, std::distance(first, last), T());
}

} // namespace multi
} // namespace boost


#endif
