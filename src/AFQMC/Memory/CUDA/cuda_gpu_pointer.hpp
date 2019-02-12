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
#include "AFQMC/Kernels/reference_operations.cuh"

#include "multi/array_ref.hpp"
#include "mpi3/shared_window.hpp"
#include "AFQMC/Memory/raw_pointers.hpp"

namespace qmc_cuda {

using qmcplusplus::afqmc::to_address;

template<class T> struct cuda_gpu_allocator;
template<class T> struct cuda_gpu_ptr;

// no const for now
template<class T>
struct cuda_gpu_reference {

  public:

  using value_type = T; 
  using decay_value_type = typename std::decay<T>::type;
  using pointer = cuda_gpu_ptr<T>;

  // must construct through a gpu_ptr for now, to keep some sort of control/safety 
  cuda_gpu_reference(pointer const& gpu_ptr) : impl_(to_address(gpu_ptr)) {}
  cuda_gpu_reference(cuda_gpu_reference<T> const& gpu_ref) = default; 
 
  // assignment
  cuda_gpu_reference& operator=(cuda_gpu_reference const& x) {
    if(cudaSuccess != cudaMemcpy(impl_,x.impl_,sizeof(T),cudaMemcpyDefault))
     throw std::runtime_error("Error: cudaMemcpy returned error code.");
    return *this;
  } 

  cuda_gpu_reference& operator=(value_type const& x) {
    if(cudaSuccess != cudaMemcpy(impl_,std::addressof(x),sizeof(T),cudaMemcpyDefault))
     throw std::runtime_error("Error: cudaMemcpy returned error code.");
    return *this;
  }

  operator value_type() const { return this->val(); }

  operator value_type&() {
    if(cudaSuccess != cudaMemcpy(std::addressof(host_impl_),impl_,sizeof(T),cudaMemcpyDefault))
     throw std::runtime_error("Error: cudaMemcpy returned error code.");
    return host_impl_;
  }

  pointer operator& () const { return pointer{impl_}; }

  void swap(cuda_gpu_reference& other) { std::swap(impl_,other.impl_); }

  cuda_gpu_reference&  operator++ (void) { kernels::op_plus(impl_,T(1)); return *this; }
  cuda_gpu_reference&  operator-- (void) { kernels::op_minus(impl_,T(1)); return *this; }
 
  value_type  operator++ (int) {
    value_type res;
    if(cudaSuccess != cudaMemcpy(std::addressof(res),impl_,sizeof(T),cudaMemcpyDefault))
     throw std::runtime_error("Error: cudaMemcpy returned error code.");
    return res++;
  }
  value_type  operator-- (int) {
    value_type res;
    if(cudaSuccess != cudaMemcpy(std::addressof(res),impl_,sizeof(T),cudaMemcpyDefault))
     throw std::runtime_error("Error: cudaMemcpy returned error code.");
    return res--;
  }

  cuda_gpu_reference&  operator+= (const value_type &rhs) 
    { kernels::op_plus(impl_,rhs); return *this; }
  cuda_gpu_reference&  operator-= (const value_type &rhs) 
    { kernels::op_minus(impl_,rhs); return *this; }
  cuda_gpu_reference&  operator*= (const value_type &rhs) 
    { kernels::op_times(impl_,rhs); return *this; }
  cuda_gpu_reference&  operator/= (const value_type &rhs) 
    { kernels::op_div(impl_,rhs); return *this; } 

  cuda_gpu_reference&  operator+= (cuda_gpu_reference const& rhs) 
    { kernels::op_plus(impl_,rhs.val()); return *this; }
  cuda_gpu_reference&  operator-= (cuda_gpu_reference const& rhs) 
    { kernels::op_minus(impl_,rhs.val()); return *this; }
  cuda_gpu_reference&  operator*= (cuda_gpu_reference const& rhs) 
    { kernels::op_times(impl_,rhs.val()); return *this; }
  cuda_gpu_reference&  operator/= (cuda_gpu_reference const& rhs) 
    { kernels::op_div(impl_,rhs.val()); return *this; }

  friend value_type& operator+=(value_type& lhs, cuda_gpu_reference const& rhs)
  { lhs += rhs.val(); return lhs; }
  friend value_type& operator-=(value_type& lhs, cuda_gpu_reference const& rhs)
  { lhs -= rhs.val(); return lhs; }
  friend value_type& operator*=(value_type& lhs, cuda_gpu_reference const& rhs)
  { lhs *= rhs.val(); return lhs; }
  friend value_type& operator/=(value_type& lhs, cuda_gpu_reference const& rhs)
  { lhs /= rhs.val(); return lhs; }

  value_type  operator+ (value_type const& rhs) const { return this->val()+rhs; }
  value_type  operator- (value_type const& rhs) const { return this->val()-rhs; }
  value_type  operator/ (value_type const& rhs) const { return this->val()/rhs; }
  value_type  operator* (value_type const& rhs) const { return this->val()*rhs; }

  value_type  operator+ (cuda_gpu_reference const& rhs) const 
    { return this->val()+rhs.val(); }
  value_type  operator- (cuda_gpu_reference const& rhs) const 
    { return this->val()-rhs.val(); }
  value_type  operator* (cuda_gpu_reference const& rhs) const 
    { return this->val()*rhs.val(); }
  value_type  operator/ (cuda_gpu_reference const& rhs) const 
    { return this->val()/rhs.val(); }

  friend value_type operator+(value_type lhs, cuda_gpu_reference const& rhs) 
  { return lhs + rhs.val(); }
  friend value_type operator-(value_type lhs, cuda_gpu_reference const& rhs) 
  { return lhs - rhs.val(); }
  friend value_type operator*(value_type lhs, cuda_gpu_reference const& rhs) 
  { return lhs * rhs.val(); }
  friend value_type operator/(value_type lhs, cuda_gpu_reference const& rhs) 
  { return lhs / rhs.val(); }

  bool operator==(value_type const& rhs) const { return this->val()==rhs; } 
  bool operator!=(value_type const& rhs) const { return this->val()!=rhs; } 
  bool operator>(value_type const& rhs) const { return this->val()>rhs; } 
  bool operator<(value_type const& rhs) const { return this->val()<rhs; } 
  bool operator>=(value_type const& rhs) const { return this->val()>=rhs; } 
  bool operator<=(value_type const& rhs) const { return this->val()<=rhs; } 

  bool operator==(cuda_gpu_reference const& rhs) const { return this->val()==rhs.val(); } 
  bool operator!=(cuda_gpu_reference const& rhs) const { return this->val()!=rhs.val(); } 
  bool operator>(cuda_gpu_reference const& rhs) const { return this->val()>rhs.val(); } 
  bool operator<(cuda_gpu_reference const& rhs) const { return this->val()<rhs.val(); } 
  bool operator>=(cuda_gpu_reference const& rhs) const { return this->val()>=rhs.val(); } 
  bool operator<=(cuda_gpu_reference const& rhs) const { return this->val()<=rhs.val(); } 

  friend bool operator==(value_type const& lhs, cuda_gpu_reference const& rhs)
    { return lhs==rhs.val(); }
  friend bool operator!=(value_type const& lhs, cuda_gpu_reference const& rhs)
    { return lhs!=rhs.val(); }
  friend bool operator>(value_type const& lhs, cuda_gpu_reference const& rhs)
    { return lhs>rhs.val(); }
  friend bool operator<(value_type const& lhs, cuda_gpu_reference const& rhs)
    { return lhs<rhs.val(); }
  friend bool operator>=(value_type const& lhs, cuda_gpu_reference const& rhs)
    { return lhs>=rhs.val(); }
  friend bool operator<=(value_type const& lhs, cuda_gpu_reference const& rhs)
    { return lhs<=rhs.val(); }

  friend std::ostream& operator<<(std::ostream& os, cuda_gpu_reference const& obj)
  {
    os<<obj.val();
    return os;
  }

  friend std::istream& operator<<(std::istream& is, cuda_gpu_reference & obj)
  {
    value_type val;
    is>>val;
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

  value_type val() const {
    decay_value_type res;
    if(cudaSuccess != cudaMemcpy(std::addressof(res),impl_,sizeof(value_type),cudaMemcpyDefault))
     throw std::runtime_error("Error: cudaMemcpy returned error code.");
    return value_type(res);
  }

};


struct base_cuda_gpu_ptr
{
  static gpu_handles handles;
};

// this class is not safe, since it allows construction of a gpu_ptr from a raw ptr
// which might not be in gpu memory. Fix this!!!
template<class T>
struct cuda_gpu_ptr: base_cuda_gpu_ptr{
  using difference_type = std::ptrdiff_t;
  using value_type = T;
  using element_type = T;
  using const_value_type = T const;
  using pointer = T*;
  using const_pointer = T const*;
// this is wrong!!! but no synthetic references yet!!!
  using reference = cuda_gpu_reference<T>;
  using const_reference = cuda_gpu_reference<T> const;
  using iterator_category = std::random_access_iterator_tag; 
  static const int memory_type = GPU_MEMORY_POINTER_TYPE; 
  using default_allocator_type = cuda_gpu_allocator<T>;
  friend class cuda_gpu_allocator<T>;
  friend class cuda_gpu_ptr<typename std::decay<T>::type>;
  default_allocator_type default_allocator() const{ return cuda_gpu_allocator<T>{}; };
  cuda_gpu_ptr() = default;
  cuda_gpu_ptr(std::nullptr_t): impl_(nullptr){}  
// eventually check if memory types and blas types are convertible, e.g. CPU_MEMORY to CPU_OUTOFCARD
  template<typename Q>
  cuda_gpu_ptr(cuda_gpu_ptr<Q> const& ptr):impl_(ptr.impl_) {}
  reference operator*() const{ return reference(cuda_gpu_ptr{impl_}); }
  reference operator[](std::ptrdiff_t n) const { return reference(cuda_gpu_ptr{impl_ + n}); }
  T* operator->() const{return impl_;}
  explicit operator bool() const{return (impl_!=nullptr);}
//  operator cuda_gpu_ptr<T const>() const{return cuda_gpu_ptr<T const>{impl_}; }
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
  cuda_gpu_ptr<Q> pointer_cast(){
    cuda_gpu_ptr<Q> res;
    res.impl_ = reinterpret_cast<Q*>(impl_);
    return res;
  }
  template<class Q>
  friend cuda_gpu_ptr<Q> pointer_cast(cuda_gpu_ptr const& self){
    return cuda_gpu_ptr<Q>{reinterpret_cast<Q*>(self.impl_)};
  }
  T* impl_;
  protected:
  cuda_gpu_ptr(T* impl__):
                         impl_(impl__) {}
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
    cuda_gpu_ptr<T> res;
    res.impl_=p;
    return res;
//    return cuda_gpu_ptr<T>{p};
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

template<typename T>
cuda_gpu_ptr<T> fill(cuda_gpu_ptr<T> first, cuda_gpu_ptr<T> last, T const& val){
  return fill_n(first,std::distance(first,last),val); 
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

template<typename T>
cuda_gpu_ptr<T> uninitialized_fill(cuda_gpu_ptr<T> first, cuda_gpu_ptr<T> last, T const& val){
  return uninitialized_fill_n(first,std::distance(first,last),val);
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

using qmcplusplus::afqmc::to_address;

/**************** copy *****************/
// Can always call cudaMemcopy2D like you do in the blas backend

template<typename T, typename Size>
multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> uninitialized_fill_n(multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> first, Size n, T const& val){
  if(n == 0) return first;
  kernels::fill_n(to_address(base(first)), n, stride(first), val);
  return first + n;
}

template<typename T, typename Size>
multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> uninitialized_fill(      
                    multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> first, 
                    multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> last, T const& val){
  assert( stride(first) == stride(last) );
  if(std::distance(first,last) == 0 ) return first;
  kernels::fill_n(to_address(base(first)), std::distance(first,last), stride(first), val);
  return first + std::distance(first,last);
}

template<class T>
multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> copy( 
           multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> first,
           multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> last,
           multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> dest ){
  assert( stride(first) == stride(last) );
  if(std::distance(first,last) == 0 ) return dest;
  if(cudaSuccess != cudaMemcpy2D(to_address(base(dest)),sizeof(T)*stride(dest),
                                 to_address(base(first)),sizeof(T)*stride(first),
                                 sizeof(T),std::distance(first,last),cudaMemcpyDeviceToDevice))
      throw std::runtime_error("Error: cudaMemcpy2D returned error code.");
  return dest+std::distance(first,last);
}

template<class T, typename Size>
multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> copy_n( 
             multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> first,
             Size N,
             multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> dest ){
  if(N==0) return dest;  
  if(cudaSuccess != cudaMemcpy2D(to_address(base(dest)),sizeof(T)*stride(dest),
                                 to_address(base(first)),sizeof(T)*stride(first),
                                 sizeof(T),N,cudaMemcpyDeviceToDevice))
      throw std::runtime_error("Error: cudaMemcpy2D returned error code.");
  return dest+N;
}

template<class T>
multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> uninitialized_copy( 
                         multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> first,
                         multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> last,
                         multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> dest ){
  assert( stride(first) == stride(last) );
  if(std::distance(first,last) == 0 ) return dest;
  if(cudaSuccess != cudaMemcpy2D(to_address(base(dest)),sizeof(T)*stride(dest),
                                 to_address(base(first)),sizeof(T)*stride(first),
                                 sizeof(T),std::distance(first,last),cudaMemcpyDeviceToDevice))
      throw std::runtime_error("Error: cudaMemcpy2D returned error code.");
  return dest+std::distance(first,last); 
}

template<class T>
multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> uninitialized_copy(
                         multi::array_iterator<T, 1, boost::mpi3::intranode::array_ptr<T>> first,
                         multi::array_iterator<T, 1, boost::mpi3::intranode::array_ptr<T>> last,
                         multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> dest ){
  assert( stride(first) == stride(last) );
  if(std::distance(first,last) == 0 ) return dest;
  if(cudaSuccess != cudaMemcpy2D(to_address(base(dest)),sizeof(T)*stride(dest),
                                 to_address(base(first)),sizeof(T)*stride(first),
                                 sizeof(T),std::distance(first,last),cudaMemcpyHostToDevice))
      throw std::runtime_error("Error: cudaMemcpy2D returned error code.");
  return dest+std::distance(first,last);
}

template<class T>
multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> uninitialized_copy(
                         multi::array_iterator<T, 1, T*> first,
                         multi::array_iterator<T, 1, T*> last,
                         multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> dest ){
  assert( stride(first) == stride(last) );
  if(std::distance(first,last) == 0 ) return dest;
  if(cudaSuccess != cudaMemcpy2D(to_address(base(dest)),sizeof(T)*stride(dest),
                                 to_address(base(first)),sizeof(T)*stride(first),
                                 sizeof(T),std::distance(first,last),cudaMemcpyHostToDevice))
      throw std::runtime_error("Error: cudaMemcpy2D returned error code.");
  return dest+std::distance(first,last);
}
template<class T>
multi::array_iterator<T, 1, T*> uninitialized_copy(
                         multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> first,
                         multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> last,
                         multi::array_iterator<T, 1, T*> dest ){
  assert( stride(first) == stride(last) );
  if(std::distance(first,last) == 0 ) return dest;
  if(cudaSuccess != cudaMemcpy2D(to_address(base(dest)),sizeof(T)*stride(dest),
                                 to_address(base(first)),sizeof(T)*stride(first),
                                 sizeof(T),std::distance(first,last),cudaMemcpyDeviceToHost))
      throw std::runtime_error("Error: cudaMemcpy2D returned error code.");
  return dest+std::distance(first,last);
}

template<class T, typename Size>
multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> uninitialized_copy_n( 
                           multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> first,
                           Size N,
                           multi::array_iterator<T, 1, qmc_cuda::cuda_gpu_ptr<T>> dest ){
  if(N==0) return dest;
  if(cudaSuccess != cudaMemcpy2D(to_address(base(dest)),sizeof(T)*stride(dest),
                                 to_address(base(first)),sizeof(T)*stride(first),
                                 sizeof(T),N,cudaMemcpyDeviceToDevice))
      throw std::runtime_error("Error: cudaMemcpy2D returned error code.");
  return dest+N;
}

}

// 
/*
namespace boost{
namespace mpi3{
namespace detail{

template<class qmc_cuda::cuda_gpu_ptr<int>> struct iterator_category{ using type = contiguous_iterator_tag; };
template<class qmc_cuda::cuda_gpu_ptr<long>> struct iterator_category{ using type = contiguous_iterator_tag; };
template<class qmc_cuda::cuda_gpu_ptr<size_t>> struct iterator_category{ using type = contiguous_iterator_tag; };
template<class qmc_cuda::cuda_gpu_ptr<float>> struct iterator_category{ using type = contiguous_iterator_tag; };
template<class qmc_cuda::cuda_gpu_ptr<double>> struct iterator_category{ using type = contiguous_iterator_tag; };
template<class qmc_cuda::cuda_gpu_ptr<std::complex<float>>> struct iterator_category{ using type = contiguous_iterator_tag; };
template<class qmc_cuda::cuda_gpu_ptr<std::complex<double>>> struct iterator_category{ using type = contiguous_iterator_tag; };

}
}
}
*/
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
