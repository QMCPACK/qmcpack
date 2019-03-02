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
#include <iostream>
#include <string>
#include <stdexcept>

#define cublasErrorCheck(ans, cause) { cublasAssert((ans), cause, __FILE__, __LINE__); }

inline void cublasAssert(cublasStatus_t code, const std::string& cause, const char *file, int line, bool abort=true)
{
  if (code != CUBLAS_STATUS_SUCCESS)
  {
    std::string cublas_error;
    switch (code)
    {
      case CUBLAS_STATUS_NOT_INITIALIZED:
        cublas_error = "CUBLAS_STATUS_NOT_INITIALIZED"; break;
      case CUBLAS_STATUS_ALLOC_FAILED:
        cublas_error = "CUBLAS_STATUS_ALLOC_FAILED"; break;
      case CUBLAS_STATUS_INVALID_VALUE:
        cublas_error = "CUBLAS_STATUS_INVALID_VALUE"; break;
      case CUBLAS_STATUS_ARCH_MISMATCH:
        cublas_error = "CUBLAS_STATUS_ARCH_MISMATCH"; break;
      case CUBLAS_STATUS_MAPPING_ERROR:
        cublas_error = "CUBLAS_STATUS_MAPPING_ERROR"; break;
      case CUBLAS_STATUS_EXECUTION_FAILED:
        cublas_error = "CUBLAS_STATUS_EXECUTION_FAILED"; break;
      case CUBLAS_STATUS_INTERNAL_ERROR:
        cublas_error = "CUBLAS_STATUS_INTERNAL_ERROR"; break;
      default:
        cublas_error = "<unknown>";
    }

    std::ostringstream err;
    err << "cublasAssert: " << cublas_error
        << ", file " << file << ", line " << line << std::endl
        << cause << std::endl;
    std::cerr << err.str();
    //if (abort) exit(code);
    throw std::runtime_error(cause);
  }
}

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
    cublasErrorCheck( cublasSgemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc), "cublasSgemm failed!" );
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
    cublasErrorCheck( cublasDgemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc), "cublasDgemm failed!");
  }
};

}
#endif // CUBLAS_CXX_H

