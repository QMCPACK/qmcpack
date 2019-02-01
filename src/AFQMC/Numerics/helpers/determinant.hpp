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
#include "AFQMC/Memory/CUDA/cuda_gpu_pointer.hpp"
#endif
#include "AFQMC/Kernels/determinant.cuh"

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

}

#if defined(QMC_CUDA)
namespace qmc_cuda{
  // using thrust for now to avoid kernels!!!
  template<class T>
  inline void determinant_from_getrf(int n, cuda_gpu_ptr<T> A, int lda, cuda_gpu_ptr<int> piv, T*res)
  {
    kernels::determinant_from_getrf_gpu(n,to_address(A),lda,to_address(piv),res);
  }

  template<class T>
  inline T determinant_from_getrf(int n, cuda_gpu_ptr<T> A, int lda, cuda_gpu_ptr<int> piv)
  {
    return kernels::determinant_from_getrf_gpu(n,to_address(A),lda,to_address(piv));
  }
}
#endif

#endif
