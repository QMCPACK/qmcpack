//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file CUDAallocator.hpp
 */
#ifndef QMCPLUSPLUS_CUDA_ALLOCATOR_H
#define QMCPLUSPLUS_CUDA_ALLOCATOR_H

#include <cstdlib>
#include <stdexcept>
#include <cuda_runtime_api.h>

namespace qmcplusplus
{
  template<typename T, size_t Align>
  struct CUDAManagedAllocator
  {
    typedef T         value_type;
    typedef size_t    size_type;
    typedef T*        pointer;
    typedef const T*  const_pointer;

    CUDAManagedAllocator() = default;
    template <class U> CUDAManagedAllocator(const CUDAManagedAllocator<U,Align>&) {}

    template <class U> struct rebind { typedef CUDAManagedAllocator<U, Align> other; };

    T* allocate(std::size_t n)
    {
      void* pt;
      cudaError_t error = cudaMallocManaged(&pt, n*sizeof(T));
      if(error!=cudaSuccess) throw std::runtime_error("Allocation failed in CUDAManagedAllocator");
      if( (size_t(pt))&(Align-1) ) throw std::runtime_error("Unaligned memory allocated in CUDAManagedAllocator");
      return static_cast<T*>(pt);
    }
    void deallocate(T* p, std::size_t)
    {
      cudaError_t error = cudaFree(p);
      if(error!=cudaSuccess) throw std::runtime_error("Deallocation failed in CUDAManagedAllocator");
    }
  };

  template <class T1, size_t Align1, class T2, size_t Align2>
  bool operator==(const CUDAManagedAllocator<T1,Align1>&, const CUDAManagedAllocator<T2,Align2>&) { return Align1==Align2; }
  template <class T1, size_t Align1, class T2, size_t Align2>
  bool operator!=(const CUDAManagedAllocator<T1,Align1>&, const CUDAManagedAllocator<T2,Align2>&) { return Align1!=Align2; }

  template<typename T>
  struct CUDAAllocator
  {
    typedef T         value_type;
    typedef size_t    size_type;
    typedef T*        pointer;
    typedef const T*  const_pointer;

    CUDAAllocator() = default;
    template <class U> CUDAAllocator(const CUDAAllocator<U>&) {}

    template <class U> struct rebind { typedef CUDAAllocator<U> other; };

    T* allocate(std::size_t n)
    {
      void* pt;
      cudaError_t error = cudaMalloc(&pt, n*sizeof(T));
      if(error!=cudaSuccess) throw std::runtime_error("Allocation failed in CUDAAllocator");
      return static_cast<T*>(pt);
    }
    void deallocate(T* p, std::size_t)
    {
      cudaError_t error = cudaFree(p);
      if(error!=cudaSuccess) throw std::runtime_error("Deallocation failed in CUDAAllocator");
    }
  };

  template <class T1, class T2>
  bool operator==(const CUDAAllocator<T1>&, const CUDAAllocator<T2>&) { return true; }
  template <class T1, class T2>
  bool operator!=(const CUDAAllocator<T1>&, const CUDAAllocator<T2>&) { return false; }
}

#endif
