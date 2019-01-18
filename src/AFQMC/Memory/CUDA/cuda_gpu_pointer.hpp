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

#include "Configuration.h"
#include<cassert>
#include <cuda_runtime.h>
#include "cublas_v2.h"
#include "cublasXt.h"
#include "AFQMC/Memory/CUDA/cuda_utilities.h"
#include "AFQMC/Kernels/fill_n.cuh"
#include "AFQMC/Kernels/uninitialized_fill_n.cuh"
#include "AFQMC/Kernels/uninitialized_copy_n.cuh"
#include "AFQMC/Kernels/copy_n_cast.cuh"
#include "AFQMC/Kernels/print.cuh"

#include "multi/array_ref.hpp"

namespace qmc_cuda {

template<class T> struct cuda_gpu_allocator;

struct base_cuda_gpu_ptr
{
  static gpu_handles handles;
};

template<class T>
struct cuda_gpu_ptr: base_cuda_gpu_ptr{
  using difference_type = std::ptrdiff_t;
  using value_type = T;
  using element_type = T;
  using const_value_type = T const;
  using pointer = T*;
  using const_pointer = T const*;
// this is wrong!!! but no synthetic references yet!!!
  using reference = T&;
  using const_reference = T const&;
  using iterator_category = std::random_access_iterator_tag; 
  static const int memory_type = GPU_MEMORY_POINTER_TYPE; 
  using default_allocator_type = cuda_gpu_allocator<T>;
  default_allocator_type default_allocator() const{ return cuda_gpu_allocator<T>{}; };
  T* impl_;
  cuda_gpu_ptr() = default;
  cuda_gpu_ptr(T* impl__):
                         impl_(impl__) {}
// eventually check if memory types and blas types are convertible, e.g. CPU_MEMORY to CPU_OUTOFCARD
  template<typename Q,
           typename = typename std::enable_if_t<cuda_gpu_ptr<Q>::memory_type == memory_type>
              >
  cuda_gpu_ptr(cuda_gpu_ptr<Q> const& ptr):impl_(ptr.impl_) {}
  T& operator*() const{return *impl_;}
  T& operator[](std::ptrdiff_t n) const{return impl_[n];}
  T* operator->() const{return impl_;}
  explicit operator bool() const{return (impl_!=nullptr);}
  operator cuda_gpu_ptr<T const>() const{return cuda_gpu_ptr<T const>{impl_}; }
  auto operator+(std::ptrdiff_t n) const{return cuda_gpu_ptr{impl_ + n};} 
  std::ptrdiff_t operator-(cuda_gpu_ptr other) const{return std::ptrdiff_t(impl_-other.impl_);}
  cuda_gpu_ptr& operator++() {++impl_; return *this;} 
  cuda_gpu_ptr& operator--() {--impl_; return *this;} 
  cuda_gpu_ptr& operator+=(std::ptrdiff_t d){impl_ += d; return *this;}
  cuda_gpu_ptr& operator-=(std::ptrdiff_t d){impl_ -= d; return *this;}
  bool operator==(cuda_gpu_ptr const& other) const{ return impl_==other.impl_; }
  bool operator!=(cuda_gpu_ptr const& other) const{ return not (*this == other); } 
  bool operator<=(cuda_gpu_ptr<T> const& other) const{
    return impl_ <= other.impl_;
  }
  T* to_address() const {return impl_;}
  friend decltype(auto) to_address(cuda_gpu_ptr const& self){return self.to_address();}
  template<class Q>
  friend cuda_gpu_ptr<Q> pointer_cast(cuda_gpu_ptr const& self){
    return cuda_gpu_ptr<Q>{reinterpret_cast<Q*>(self.impl_)};
  }
};

//static size_t TotalGPUAlloc=0;
  
/*
 * Incomplete!!! Need to fix construct and destroy
 */
template<class T> struct cuda_gpu_allocator{
  template<class U> struct rebind{typedef cuda_gpu_allocator<U> other;};
  using element_type = T;
  using value_type = T;
  using const_value_type = T const;
  using pointer = cuda_gpu_ptr<T>;
  using const_pointer = cuda_gpu_ptr<T const>;
  using reference = T&;
  using const_reference = T const&;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

  cuda_gpu_allocator() = default; 
  ~cuda_gpu_allocator() = default;
  cuda_gpu_allocator(cuda_gpu_allocator const& other) = default;
  template<class U>
  cuda_gpu_allocator(cuda_gpu_allocator<U> const& other) {}

  cuda_gpu_ptr<T> allocate(size_type n, const void* hint = 0){
    if(n == 0) return cuda_gpu_ptr<T>{};
    T* p;
    if(cudaSuccess != cudaMalloc ((void**)&p,n*sizeof(T))) {
      std::cerr<<" Error allocating " <<n*sizeof(T)/1024.0/1024.0 <<" MBs on GPU." <<std::endl;
      throw std::runtime_error("Error: cudaMalloc returned error code."); 
    }
    return cuda_gpu_ptr<T>{p};
  }
  void deallocate(cuda_gpu_ptr<T> ptr, size_type){
    cudaFree(ptr.impl_); 
  }
  bool operator==(cuda_gpu_allocator const& other) const{
    return true; 
  }
  bool operator!=(cuda_gpu_allocator const& other) const{
    return false; 
  }
  template<class U, class... Args>
  void construct(U* p, Args&&... args){
    //::new((void*)p) U(std::forward<Args>(args)...);
  }
  template< class U >
  void destroy(U* p){
    //p->~U();
  }
};


/* Don't know how to implement this on the kernel side, without propagating the 
 * cuda code upstream due to the template needed to pass the UnaryOperator
template<class T, class F>
F for_each(cuda_gpu_ptr<T> first, cuda_gpu_ptr<T> last, F f){
        if(first == last) return f;
        return kernels::for_each(to_address(first), to_address(last), f);
}
*/
/**************** copy_n *****************/
template<typename T, typename Size>
cuda_gpu_ptr<T> copy_n(cuda_gpu_ptr<T> const A, Size n, cuda_gpu_ptr<T> B) {
  if(cudaSuccess != cudaMemcpy(to_address(B),to_address(A),n*sizeof(T),cudaMemcpyDefault))
   throw std::runtime_error("Error: cudaMemcpy returned error code.");
  return B+n;
}

template<typename T, typename Size>
cuda_gpu_ptr<T> copy_n(T* const A, Size n, cuda_gpu_ptr<T> B) {
  if(cudaSuccess != cudaMemcpy(to_address(B),A,n*sizeof(T),cudaMemcpyDefault))
   throw std::runtime_error("Error: cudaMemcpy returned error code.");
  return B+n;
}

template<typename T, typename Size>
T* copy_n(cuda_gpu_ptr<T> const A, Size n, T* B) {
  if(cudaSuccess != cudaMemcpy(B,to_address(A),n*sizeof(T),cudaMemcpyDefault))
   throw std::runtime_error("Error: cudaMemcpy returned error code.");
  return B+n;
}

/**************** copy *****************/
template<typename T>
cuda_gpu_ptr<T> copy(cuda_gpu_ptr<T> const Abeg, cuda_gpu_ptr<T> const Aend, cuda_gpu_ptr<T> B) {
  return copy_n(Abeg,std::distance(Abeg,Aend),B);
}

template<typename T>
cuda_gpu_ptr<T> copy(T* const Abeg, T* const Aend, cuda_gpu_ptr<T> B) {
  return copy_n(Abeg,std::distance(Abeg,Aend),B);
}

template<typename T>
T* copy(cuda_gpu_ptr<T> const Abeg, cuda_gpu_ptr<T> const Aend, T* B) {
  return copy_n(Abeg,std::distance(Abeg,Aend),B);
}

/**************** copy_n_cast *****************/
template<typename T, typename Q, typename Size>
cuda_gpu_ptr<T> copy_n_cast(cuda_gpu_ptr<T> const A, Size n, cuda_gpu_ptr<Q> B) {
  kernels::copy_n_cast(to_address(A),n,to_address(B));
  return B+n;
}

template<typename T, typename Size>
cuda_gpu_ptr<T> copy_n_cast(T* const A, Size n, cuda_gpu_ptr<T> B) {
  throw std::runtime_error(" Error: copy_n_cast(gpu_ptr,n,T*) is disabled.");
  return B+n;
}

template<typename T, typename Size>
T* copy_n_cast(cuda_gpu_ptr<T> const A, Size n, T* B) {
  throw std::runtime_error(" Error: copy_n_cast(gpu_ptr,n,T*) is disabled.");
  return B+n;
}

/**************** fill_n *****************/
//template<typename T, typename Size, typename... Args>
//cuda_gpu_ptr<T> fill_n(cuda_gpu_ptr<T> first, Size n, Args&&...args){
//  if(n == 0) return first;
//  kernels::fill_n(to_address(first), n, std::forward<Args>(args)...);
//  return first + n;
//}

template<typename T, typename Size>
cuda_gpu_ptr<T> fill_n(cuda_gpu_ptr<T> first, Size n, T const& val){
  if(n == 0) return first;
  kernels::fill_n(to_address(first), n, val);
  return first + n;
}

/**************** uninitialized_fill_n *****************/
/*
template<typename T, typename Size, typename... Args>
cuda_gpu_ptr<T> uninitialized_fill_n(cuda_gpu_ptr<T> first, Size n, Args&&...args){
  if(n == 0) return first;
  kernels::uninitialized_fill_n(to_address(first), n, std::forward<Args>(args)...); 
  return first + n;
}
*/

template<typename T, typename Size>
cuda_gpu_ptr<T> uninitialized_fill_n(cuda_gpu_ptr<T> first, Size n, T const& val){
  if(n == 0) return first;
  kernels::uninitialized_fill_n(to_address(first), n, val);
  return first + n;
}


/******************/


template<typename T, typename Size>
cuda_gpu_ptr<T> uninitialized_default_construct_n(cuda_gpu_ptr<T> first, Size n){
  return first+n;
// what to do???
}
template<typename T, typename Size>
cuda_gpu_ptr<T> uninitialized_value_construct_n(cuda_gpu_ptr<T> first, Size n){
  return first+n;
// what to do???
}

/**************** uninitialized_copy_n *****************/
template<typename T, typename Size> 
cuda_gpu_ptr<T> uninitialized_copy_n(cuda_gpu_ptr<T> first, Size n, cuda_gpu_ptr<T> dest){
  if(n == 0) return dest;
  kernels::uninitialized_copy_n(n,to_address(first), 1, to_address(dest), 1);
  return dest + n;
}

template<class T> 
cuda_gpu_ptr<T> uninitialized_copy(cuda_gpu_ptr<T> first, cuda_gpu_ptr<T> last, cuda_gpu_ptr<T> dest){
  return uninitialized_copy_n(first,std::distance(first,last),dest); 
}


/**************** destroy_n *****************/
// NOTE: Not sure what to do here
// should at least guard agains non-trivial types
template<typename T, typename Size>
cuda_gpu_ptr<T> destroy_n(cuda_gpu_ptr<T> first, Size n){
  return first + n;
}

/**************** print *****************/
template<typename T>
void print(std::string str, cuda_gpu_ptr<T> p, int n) {
  kernels::print(str,to_address(p),n);
}

}

namespace boost::multi{

/**************** copy *****************/
// Can always call cudaMemcopy2D like you do in the blas backend
template<class T>
multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> copy( 
           multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> first,
           multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> last,
           multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> dest ){
  assert( stride(first) == stride(last) );
  if(std::distance(first,last) == 0 ) return dest;
  kernels::uninitialized_copy_n(std::distance(first,last),to_address(base(first)),stride(first),
                                                          to_address(base(dest)),stride(dest));
  return dest+std::distance(first,last);
}

template<class T, typename Size>
multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> copy_n( 
             multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> first,
             Size N,
             multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> dest ){
  if(N==0) return dest;  
  kernels::uninitialized_copy_n(N,to_address(base(first)),stride(first),
                                  to_address(base(dest)),stride(dest));
  return dest+N;
}

template<class T>
multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> uninitialized_copy( 
                         multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> first,
                         multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> last,
                         multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> dest ){
  assert( stride(first) == stride(last) );
  if(std::distance(first,last) == 0 ) return dest;
  kernels::uninitialized_copy_n(std::distance(first,last),to_address(base(first)),stride(first),
                                                          to_address(base(dest)),stride(dest));
  return dest+std::distance(first,last); 
}

template<class T, typename Size>
multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> uninitialized_copy_n( 
                           multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> first,
                           Size N,
                           multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> dest ){
  if(N==0) return dest;
  kernels::uninitialized_copy_n(N,to_address(base(first)),stride(first),
                                  to_address(base(dest)),stride(dest));
  return dest+N;
}

}

/*
namespace boost
{
namespace multi
{
  template<>
  struct pointer_traits<qmc_cuda::cuda_gpu_ptr<std::complex<double>>>:  
            std::pointer_traits<std::complex<double>>{
//    using element_type = std::complex<double>;
    using allocator_type = qmc_cuda::cuda_gpu_allocator<std::complex<double>>;
    static allocator_type allocator_of(qmc_cuda::cuda_gpu_ptr<std::complex<double>>){return allocator_type{};}
  };
  template<>
  struct pointer_traits<qmc_cuda::cuda_gpu_ptr<std::complex<float>>> {
    using element_type = std::complex<float>;
    using allocator_type = qmc_cuda::cuda_gpu_allocator<std::complex<float>>;
    static allocator_type allocator_of(qmc_cuda::cuda_gpu_ptr<std::complex<float>>){return allocator_type{};}
  };
  template<>
  struct pointer_traits<qmc_cuda::cuda_gpu_ptr<double>> {
    using element_type = double; 
    using allocator_type = qmc_cuda::cuda_gpu_allocator<double>;
    static allocator_type allocator_of(qmc_cuda::cuda_gpu_ptr<double>){return allocator_type{};}
  };
  template<>
  struct pointer_traits<qmc_cuda::cuda_gpu_ptr<float>> {
    using element_type = float; 
    using allocator_type = qmc_cuda::cuda_gpu_allocator<float>;
    static allocator_type allocator_of(qmc_cuda::cuda_gpu_ptr<float>){return allocator_type{};}
  };
  template<>
  struct pointer_traits<qmc_cuda::cuda_gpu_ptr<int>> {
    using element_type = int; 
    using allocator_type = qmc_cuda::cuda_gpu_allocator<int>;
    static allocator_type allocator_of(qmc_cuda::cuda_gpu_ptr<int>){return allocator_type{};}
  };
}
}
*/
  
#endif
