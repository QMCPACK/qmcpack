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

#ifndef QMCPLUSPLUS_CUBLAS_H
#define QMCPLUSPLUS_CUBLAS_H

#include <cublas_v2.h>
#include <complex>

namespace qmcplusplus {

struct cuBLAS
{
  static inline
  void gemm(cublasHandle_t& handle,
            const cublasOperation_t& transa, const cublasOperation_t& transb,
            int m, int n, int k,
            const float *alpha,
            const float *A, int lda,
            const float *B, int ldb,
            const float *beta,
            float *C, int ldc)
  {
    cublasStatus_t cublas_error;
    cublas_error = cublasSgemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    if( cublas_error != CUBLAS_STATUS_SUCCESS ) throw std::runtime_error("cublasSgemm failed error = " + cublas_error);
  }

  static inline
  void gemm(cublasHandle_t& handle,
            const cublasOperation_t& transa, const cublasOperation_t& transb,
            int m, int n, int k,
            const double *alpha,
            const double *A, int lda,
            const double *B, int ldb,
            const double *beta,
            double *C, int ldc)
  {
    cublasStatus_t cublas_error;
    cublas_error = cublasDgemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    if( cublas_error != CUBLAS_STATUS_SUCCESS ) throw std::runtime_error("cublasDgemm failed error = " + cublas_error);
  }
};

}
#endif // CUBLAS_CXX_H

