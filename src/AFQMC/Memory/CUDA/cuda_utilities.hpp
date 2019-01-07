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

#ifndef AFQMC_CUDA_UTILITIES_HPP 
#define AFQMC_CUDA_UTILITIES_HPP

#define GPU_MEMORY_POINTER_TYPE      1001
#define MANAGED_MEMORY_POINTER_TYPE  2001
#define CPU_OUTOFCARS_POINTER_TYPE   3001

#include<cassert>
#include<cstdlib>
#include<iostream>
#include<stdexcept>
#include <cuda_runtime.h>
#include "cublas_v2.h"
#include "cublasXt.h"
#include "cusolverDn.h"
#include "curand.h"

namespace qmc_cuda {

  inline void cuda_check(cudaError_t sucess, std::string message="")
  {
    if(cudaSuccess != sucess) {
      std::cerr<<message <<std::endl;
      std::cerr.flush();
      throw std::runtime_error(" Error code returned by cuda. \n");
    }
  }

  inline void cublas_check(cublasStatus_t sucess, std::string message="")
  {
    if(CUBLAS_STATUS_SUCCESS != sucess) {
      std::cerr<<message <<std::endl;
      std::cerr.flush();
      throw std::runtime_error(" Error code returned by cublas. \n");
    }
  }

  inline void curand_check(curandStatus_t sucess, std::string message="")
  {
    if(CURAND_STATUS_SUCCESS != sucess) {
      std::cerr<<message <<std::endl;
      std::cerr.flush();
      throw std::runtime_error(" Error code returned by curand. \n");
    }
  }

  inline void cusolver_check(cusolverStatus_t sucess, std::string message="")
  {
    if(CUSOLVER_STATUS_SUCCESS != sucess) {
      std::cerr<<message <<std::endl;
      std::cerr.flush();
      throw std::runtime_error(" Error code returned by cusolver. \n");
    }
  }

  inline cublasOperation_t cublasOperation(char A) {
    if(A=='N' or A=='n')
      return CUBLAS_OP_N;
    else if(A=='T' or A=='t')
      return CUBLAS_OP_T;
    else if(A=='C' or A=='c')
      return CUBLAS_OP_C;
    else
      throw std::runtime_error("unknown cublasOperation option"); 
    return CUBLAS_OP_N;
  }

  struct gpu_handles {
    cublasHandle_t* cublas_handle;
    cublasXtHandle_t* cublasXt_handle;
    cusolverDnHandle_t* cusolverDn_handle; 
  };
}

#endif
