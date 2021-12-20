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

#ifndef CUBLASXT_FUNCTIONDEFS_H
#define CUBLASXT_FUNCTIONDEFS_H

#include <cassert>
#include <cuda_runtime.h>
#include "cublas_v2.h"
#include "cublasXt.h"
#include "AFQMC/Memory/CUDA/cuda_utilities.h"

namespace cublas
{
using qmc_cuda::cublasOperation;

// cublasXt Level 3
inline cublasStatus_t cublasXt_gemm(cublasXtHandle_t handle,
                                    char Atrans,
                                    char Btrans,
                                    int M,
                                    int N,
                                    int K,
                                    const float alpha,
                                    const float* A,
                                    int lda,
                                    const float* B,
                                    int ldb,
                                    const float beta,
                                    float* C,
                                    int ldc)
{
  cublasStatus_t success = cublasXtSgemm(handle, cublasOperation(Atrans), cublasOperation(Btrans), M, N, K, &alpha, A,
                                        lda, B, ldb, &beta, C, ldc);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublasXt_gemm(cublasXtHandle_t handle,
                                    char Atrans,
                                    char Btrans,
                                    int M,
                                    int N,
                                    int K,
                                    const double alpha,
                                    const double* A,
                                    int lda,
                                    const double* B,
                                    int ldb,
                                    const double beta,
                                    double* C,
                                    int ldc)
{
  cublasStatus_t success = cublasXtDgemm(handle, cublasOperation(Atrans), cublasOperation(Btrans), M, N, K, &alpha, A,
                                        lda, B, ldb, &beta, C, ldc);
  /*
std::cout<<" Dgemm error message " <<success <<std::endl;
using std::cout;
using std::endl;
switch(success)
{
  case CUBLAS_STATUS_NOT_INITIALIZED:
    std::cout<<"CUBLAS_STATUS_NOT_INITIALIZED";
    break;
  case CUBLAS_STATUS_ALLOC_FAILED:
    cout<<"CUBLAS_STATUS_ALLOC_FAILED";
    break;
  case CUBLAS_STATUS_INVALID_VALUE:
    cout<<"CUBLAS_STATUS_INVALID_VALUE";
    break;
  case CUBLAS_STATUS_EXECUTION_FAILED:
    cout<<"CUBLAS_STATUS_EXECUTION_FAILED";
    break;
}
std::cout<<std::endl;
*/
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublasXt_gemm(cublasXtHandle_t handle,
                                    char Atrans,
                                    char Btrans,
                                    int M,
                                    int N,
                                    int K,
                                    const std::complex<float> alpha,
                                    const std::complex<float>* A,
                                    int lda,
                                    const std::complex<float>* B,
                                    int ldb,
                                    const std::complex<float> beta,
                                    std::complex<float>* C,
                                    int ldc)
{
  cublasStatus_t success =
      cublasXtCgemm(handle, cublasOperation(Atrans), cublasOperation(Btrans), M, N, K,
                    reinterpret_cast<cuComplex const*>(&alpha), reinterpret_cast<cuComplex const*>(A), lda,
                    reinterpret_cast<cuComplex const*>(B), ldb, reinterpret_cast<cuComplex const*>(&beta),
                    reinterpret_cast<cuComplex*>(C), ldc);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublasXt_gemm(cublasXtHandle_t handle,
                                    char Atrans,
                                    char Btrans,
                                    int M,
                                    int N,
                                    int K,
                                    const std::complex<double> alpha,
                                    const std::complex<double>* A,
                                    int lda,
                                    const std::complex<double>* B,
                                    int ldb,
                                    const std::complex<double> beta,
                                    std::complex<double>* C,
                                    int ldc)
{
  cublasStatus_t success =
      cublasXtZgemm(handle, cublasOperation(Atrans), cublasOperation(Btrans), M, N, K,
                    reinterpret_cast<cuDoubleComplex const*>(&alpha), reinterpret_cast<cuDoubleComplex const*>(A), lda,
                    reinterpret_cast<cuDoubleComplex const*>(B), ldb, reinterpret_cast<cuDoubleComplex const*>(&beta),
                    reinterpret_cast<cuDoubleComplex*>(C), ldc);
  cudaDeviceSynchronize();
  return success;
}

} // namespace cublas

#endif
