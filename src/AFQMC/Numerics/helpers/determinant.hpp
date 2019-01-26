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

#ifndef AFQMC_NUMERICS_HELPERS_HPP
#define AFQMC_NUMERICS_HELPERS_HPP

#include<cassert>
#if defined(QMC_CUDA)
#include <boost/stacktrace.hpp>
#include "AFQMC/Memory/CUDA/cuda_gpu_pointer.hpp"
#include "AFQMC/Kernels/determinant.cuh"
#endif

namespace ma 
{
  template<class T>
  inline T determinant_from_getrf(int n, T* M, int lda, int* pivot)
  {
    T res(1.0);
    for(int i=0, ip=1; i != n; i++, ip++){
      if(pivot[i]==ip){
        res *= +static_cast<T>(M[i*lda+i]);
      }else{
        res *= -static_cast<T>(M[i*lda+i]);
      }
    }
    return res;
  } 

  template<class T>
  inline void determinant_from_getrf(int n, T* M, int lda, int* pivot, T* res)
  {
    *res = T(1.0);
    for(int i=0, ip=1; i != n; i++, ip++){
      if(pivot[i]==ip){
        *res *= +static_cast<T>(M[i*lda+i]);
      }else{
        *res *= -static_cast<T>(M[i*lda+i]);
      }
    }
  }

  template<class T>
  inline void determinant_from_geqrf(int n, T* M, int lda, T* buff, T* res)
  {
    *res = T(1.0); 
    for (int i = 0; i < n; i++) { 
      if (real(M[i*lda+i]) < 0) 
        buff[i]=T(-1.0); 
      else 
        buff[i]=T(1.0); 
      *res *= buff[i]*M[i*lda+i];
    }
  }

  template<class T>
  inline void scale_columns(int n, int m, T* A, int lda, T* scl)
  {
    for (int i = 0; i < n; i++) 
      for (int j = 0; j < m; j++) 
        A[ i*lda + j ] *= scl[j];
  }

}

#if defined(QMC_CUDA)
namespace qmc_cuda{
  // using thrust for now to avoid kernels!!!
  template<class T>
  inline void determinant_from_getrf(int n, cuda_gpu_ptr<T> A, int lda, cuda_gpu_ptr<int> piv, T* res)
  {
    kernels::determinant_from_getrf_gpu(n,to_address(A),lda,to_address(piv),res);
  }

  template<class T>
  inline T determinant_from_getrf(int n, cuda_gpu_ptr<T> A, int lda, cuda_gpu_ptr<int> piv)
  {
    return kernels::determinant_from_getrf_gpu(n,to_address(A),lda,to_address(piv));
  }

  template<class T>
  inline void determinant_from_geqrf(int n, cuda_gpu_ptr<T> M, int lda, cuda_gpu_ptr<T> buff, T* res)
  {
    return kernels::determinant_from_geqrf_gpu(n,to_address(M),lda,to_address(buff),res);
  }

  template<class T>
  inline void scale_columns(int n, int m, cuda_gpu_ptr<T> A, int lda, cuda_gpu_ptr<T> scl)
  {
    return kernels::scale_columns(n,m,to_address(A),lda,to_address(scl));
  }

  template<class ptrA, class ptrB>
  inline void scale_columns(int n, int m, ptrA A, int lda, ptrB scl)
  {
    std::cout << boost::stacktrace::stacktrace();
    throw std::runtime_error("Error: Calling qmc_cuda::scale_columns atch all.");
  }

}
#endif

#endif
