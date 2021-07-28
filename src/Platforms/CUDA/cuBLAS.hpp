//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
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
#include "CUDATypeMapping.hpp"
#include "type_traits/type_manipulation.hpp"

#define cublasErrorCheck(ans, cause)                \
  {                                                 \
    cublasAssert((ans), cause, __FILE__, __LINE__); \
  }

/// prints cuBLAS error messages. Always use cublasErrorCheck macro.
inline void cublasAssert(cublasStatus_t code, const std::string& cause, const char* file, int line, bool abort = true)
{
  if (code != CUBLAS_STATUS_SUCCESS)
  {
    std::string cublas_error;
    switch (code)
    {
    case CUBLAS_STATUS_NOT_INITIALIZED:
      cublas_error = "CUBLAS_STATUS_NOT_INITIALIZED";
      break;
    case CUBLAS_STATUS_ALLOC_FAILED:
      cublas_error = "CUBLAS_STATUS_ALLOC_FAILED";
      break;
    case CUBLAS_STATUS_INVALID_VALUE:
      cublas_error = "CUBLAS_STATUS_INVALID_VALUE";
      break;
    case CUBLAS_STATUS_ARCH_MISMATCH:
      cublas_error = "CUBLAS_STATUS_ARCH_MISMATCH";
      break;
    case CUBLAS_STATUS_MAPPING_ERROR:
      cublas_error = "CUBLAS_STATUS_MAPPING_ERROR";
      break;
    case CUBLAS_STATUS_EXECUTION_FAILED:
      cublas_error = "CUBLAS_STATUS_EXECUTION_FAILED";
      break;
    case CUBLAS_STATUS_INTERNAL_ERROR:
      cublas_error = "CUBLAS_STATUS_INTERNAL_ERROR";
      break;
    case CUBLAS_STATUS_NOT_SUPPORTED:
      cublas_error = "CUBLAS_STATUS_NOT_SUPPORTED";
      break;
    case CUBLAS_STATUS_LICENSE_ERROR:
      cublas_error = "CUBLAS_STATUS_LICENSE_ERROR";
      break;
    default:
      cublas_error = "<unknown>";
    }

    std::ostringstream err;
    err << "cublasAssert: " << cublas_error << ", file " << file << " , line " << line << std::endl
        << cause << std::endl;
    std::cerr << err.str();
    //if (abort) exit(code);
    throw std::runtime_error(cause);
  }
}

namespace qmcplusplus
{
/** interface to cuBLAS calls for different data types S/C/D/Z
 */
namespace cuBLAS
{
inline cublasStatus_t copy(cublasHandle_t& handle,
                           int n,
                           const float* x,
                           int incx,
                           float* y,
                           int incy)
{
  return cublasScopy(handle, n, x, incx, y, incy);
}

inline cublasStatus_t copy(cublasHandle_t& handle,
                           int n,
                           const double* x,
                           int incx,
                           double* y,
                           int incy)
{
  return cublasDcopy(handle, n, x, incx, y, incy);
}

inline cublasStatus_t copy(cublasHandle_t& handle,
                           int n,
                           const std::complex<float>* x,
                           int incx,
                           std::complex<float>* y,
                           int incy)
{
  return cublasCcopy(handle, n, castCUDAType(x), incx, castCUDAType(y), incy);
}

inline cublasStatus_t copy(cublasHandle_t& handle,
                           int n,
                           const std::complex<double>* x,
                           int incx,
                           std::complex<double>* y,
                           int incy)
{
  return cublasZcopy(handle, n, castCUDAType(x), incx, castCUDAType(y), incy);
}

inline cublasStatus_t geam(cublasHandle_t& handle,
                           cublasOperation_t transa,
                           cublasOperation_t transb,
                           int m,
                           int n,
                           const float* alpha,
                           const float* A,
                           int lda,
                           const float* beta,
                           const float* B,
                           int ldb,
                           float* C,
                           int ldc)
{
  return cublasSgeam(handle, transa, transb, m, n, alpha, A, lda, beta, B, ldb, C, ldc);
}

inline cublasStatus_t geam(cublasHandle_t& handle,
                           cublasOperation_t transa,
                           cublasOperation_t transb,
                           int m,
                           int n,
                           const double* alpha,
                           const double* A,
                           int lda,
                           const double* beta,
                           const double* B,
                           int ldb,
                           double* C,
                           int ldc)
{
  return cublasDgeam(handle, transa, transb, m, n, alpha, A, lda, beta, B, ldb, C, ldc);
}

inline cublasStatus_t geam(cublasHandle_t& handle,
                           cublasOperation_t transa,
                           cublasOperation_t transb,
                           int m,
                           int n,
                           const std::complex<double>* alpha,
                           const std::complex<double>* A,
                           int lda,
                           const std::complex<double>* beta,
                           const std::complex<double>* B,
                           int ldb,
                           std::complex<double>* C,
                           int ldc)
{
  return cublasZgeam(handle, transa, transb, m, n, castCUDAType(alpha), castCUDAType(A), lda, castCUDAType(beta),
                     castCUDAType(B), ldb, castCUDAType(C), ldc);
}

inline cublasStatus_t geam(cublasHandle_t& handle,
                           cublasOperation_t transa,
                           cublasOperation_t transb,
                           int m,
                           int n,
                           const std::complex<float>* alpha,
                           const std::complex<float>* A,
                           int lda,
                           const std::complex<float>* beta,
                           const std::complex<float>* B,
                           int ldb,
                           std::complex<float>* C,
                           int ldc)
{
  return cublasCgeam(handle, transa, transb, m, n, castCUDAType(alpha), castCUDAType(A), lda, castCUDAType(beta),
                     castCUDAType(B), ldb, castCUDAType(C), ldc);
}

inline cublasStatus_t gemm(cublasHandle_t& handle,
                           const cublasOperation_t transa,
                           const cublasOperation_t transb,
                           int m,
                           int n,
                           int k,
                           const float* alpha,
                           const float* A,
                           int lda,
                           const float* B,
                           int ldb,
                           const float* beta,
                           float* C,
                           int ldc)
{
  return cublasSgemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

inline cublasStatus_t gemm(cublasHandle_t& handle,
                           const cublasOperation_t& transa,
                           const cublasOperation_t& transb,
                           int m,
                           int n,
                           int k,
                           const std::complex<float>* alpha,
                           const std::complex<float>* A,
                           int lda,
                           const std::complex<float>* B,
                           int ldb,
                           const std::complex<float>* beta,
                           std::complex<float>* C,
                           int ldc)
{
  return cublasCgemm(handle, transa, transb, m, n, k, castCUDAType(alpha), castCUDAType(A), lda, castCUDAType(B), ldb,
                     castCUDAType(beta), castCUDAType(C), ldc);
}

inline cublasStatus_t gemm(cublasHandle_t& handle,
                           const cublasOperation_t& transa,
                           const cublasOperation_t& transb,
                           int m,
                           int n,
                           int k,
                           const double* alpha,
                           const double* A,
                           int lda,
                           const double* B,
                           int ldb,
                           const double* beta,
                           double* C,
                           int ldc)
{
  return cublasDgemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

inline cublasStatus_t gemm(cublasHandle_t& handle,
                           const cublasOperation_t& transa,
                           const cublasOperation_t& transb,
                           int m,
                           int n,
                           int k,
                           const std::complex<double>* alpha,
                           const std::complex<double>* A,
                           int lda,
                           const std::complex<double>* B,
                           int ldb,
                           const std::complex<double>* beta,
                           std::complex<double>* C,
                           int ldc)
{
  return cublasZgemm(handle, transa, transb, m, n, k, castCUDAType(alpha), castCUDAType(A), lda, castCUDAType(B), ldb,
                     castCUDAType(beta), castCUDAType(C), ldc);
}

inline cublasStatus_t gemv(cublasHandle_t& handle,
                           const cublasOperation_t& transa,
                           int m,
                           int n,
                           const float* alpha,
                           const float* A,
                           int lda,
                           const float* x,
                           int incx,
                           const float* beta,
                           float* y,
                           int incy)
{
  return cublasSgemv(handle, transa, m, n, alpha, A, lda, x, incx, beta, y, incy);
}

inline cublasStatus_t gemv(cublasHandle_t& handle,
                           const cublasOperation_t& transa,
                           int m,
                           int n,
                           const double* alpha,
                           const double* A,
                           int lda,
                           const double* x,
                           int incx,
                           const double* beta,
                           double* y,
                           int incy)
{
  return cublasDgemv(handle, transa, m, n, alpha, A, lda, x, incx, beta, y, incy);
}

inline cublasStatus_t gemv(cublasHandle_t& handle,
                           const cublasOperation_t& transa,
                           int m,
                           int n,
                           const std::complex<float>* alpha,
                           const std::complex<float>* A,
                           int lda,
                           const std::complex<float>* x,
                           int incx,
                           const std::complex<float>* beta,
                           std::complex<float>* y,
                           int incy)
{
  return cublasCgemv(handle, transa, m, n, castCUDAType(alpha), castCUDAType(A), lda, castCUDAType(x), incx,
                     castCUDAType(beta), castCUDAType(y), incy);
}

inline cublasStatus_t gemv(cublasHandle_t& handle,
                           const cublasOperation_t& transa,
                           int m,
                           int n,
                           const std::complex<double>* alpha,
                           const std::complex<double>* A,
                           int lda,
                           const std::complex<double>* x,
                           int incx,
                           const std::complex<double>* beta,
                           std::complex<double>* y,
                           int incy)
{
  return cublasZgemv(handle, transa, m, n, castCUDAType(alpha), castCUDAType(A), lda, castCUDAType(x), incx,
                     castCUDAType(beta), castCUDAType(y), incy);
}

inline cublasStatus_t gemm_batched(cublasHandle_t& handle,
                                   const cublasOperation_t& transa,
                                   const cublasOperation_t& transb,
                                   int m,
                                   int n,
                                   int k,
                                   const float* alpha,
                                   const float* const A[],
                                   int lda,
                                   const float* const B[],
                                   int ldb,
                                   const float* beta,
                                   float* const C[],
                                   int ldc,
                                   int batchCount)
{
  return cublasSgemmBatched(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, batchCount);
}

inline cublasStatus_t gemm_batched(cublasHandle_t& handle,
                                   const cublasOperation_t& transa,
                                   const cublasOperation_t& transb,
                                   int m,
                                   int n,
                                   int k,
                                   const std::complex<float>* alpha,
                                   const std::complex<float>* const A[],
                                   int lda,
                                   const std::complex<float>* const B[],
                                   int ldb,
                                   const std::complex<float>* beta,
                                   std::complex<float>* const C[],
                                   int ldc,
                                   int batchCount)
{
  // This is necessary to not break the complex CUDA type mapping semantics while
  // dealing with the const cuComplex * A[] style API of cuBLAS
  // C++ makes you jump through some hoops to remove the bottom const on a double pointer.
  // see typetraits/type_manipulation.hpp
  auto non_const_A = const_cast<BottomConstRemoved<decltype(A)>::type>(A);
  auto non_const_B = const_cast<BottomConstRemoved<decltype(B)>::type>(B);
  auto non_const_C = const_cast<BottomConstRemoved<decltype(C)>::type>(C);

  return cublasCgemmBatched(handle, transa, transb, m, n, k, castCUDAType(alpha), castCUDAType(non_const_A), lda,
                            castCUDAType(non_const_B), ldb, castCUDAType(beta), castCUDAType(non_const_C), ldc,
                            batchCount);
}

inline cublasStatus_t gemm_batched(cublasHandle_t& handle,
                                   const cublasOperation_t& transa,
                                   const cublasOperation_t& transb,
                                   int m,
                                   int n,
                                   int k,
                                   const double* alpha,
                                   const double* const A[],
                                   int lda,
                                   const double* const B[],
                                   int ldb,
                                   const double* beta,
                                   double* const C[],
                                   int ldc,
                                   int batchCount)
{
  return cublasDgemmBatched(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, batchCount);
}

inline cublasStatus_t gemm_batched(cublasHandle_t& handle,
                                   const cublasOperation_t& transa,
                                   const cublasOperation_t& transb,
                                   int m,
                                   int n,
                                   int k,
                                   const std::complex<double>* alpha,
                                   const std::complex<double>* const A[],
                                   int lda,
                                   const std::complex<double>* const B[],
                                   int ldb,
                                   const std::complex<double>* beta,
                                   std::complex<double>* const C[],
                                   int ldc,
                                   int batchCount)
{
  auto non_const_A = const_cast<BottomConstRemoved<decltype(A)>::type>(A);
  auto non_const_B = const_cast<BottomConstRemoved<decltype(B)>::type>(B);
  auto non_const_C = const_cast<BottomConstRemoved<decltype(C)>::type>(C);

  return cublasZgemmBatched(handle, transa, transb, m, n, k, castCUDAType(alpha), castCUDAType(non_const_A), lda,
                            castCUDAType(non_const_B), ldb, castCUDAType(beta), castCUDAType(non_const_C), ldc,
                            batchCount);
}

inline cublasStatus_t ger(cublasHandle_t& handle,
                          int m,
                          int n,
                          const float* alpha,
                          const float* x,
                          int incx,
                          const float* y,
                          int incy,
                          float* const A,
                          int lda)
{
  return cublasSger(handle, m, n, alpha, x, incx, y, incy, A, lda);
}

inline cublasStatus_t ger(cublasHandle_t& handle,
                          int m,
                          int n,
                          const double* alpha,
                          const double* x,
                          int incx,
                          const double* y,
                          int incy,
                          double* const A,
                          int lda)
{
  return cublasDger(handle, m, n, alpha, x, incx, y, incy, A, lda);
}

inline cublasStatus_t ger(cublasHandle_t& handle,
                          int m,
                          int n,
                          const std::complex<float>* alpha,
                          const std::complex<float>* x,
                          int incx,
                          const std::complex<float>* y,
                          int incy,
                          std::complex<float>* A,
                          int lda)
{
  return cublasCgerc(handle, m, n, castCUDAType(alpha), castCUDAType(x), incx, castCUDAType(y), incy, castCUDAType(A),
                     lda);
}

inline cublasStatus_t ger(cublasHandle_t& handle,
                          int m,
                          int n,
                          const std::complex<double>* alpha,
                          const std::complex<double>* x,
                          int incx,
                          const std::complex<double>* y,
                          int incy,
                          std::complex<double>* A,
                          int lda)
{
  return cublasZgerc(handle, m, n, castCUDAType(alpha), castCUDAType(x), incx, castCUDAType(y), incy, castCUDAType(A),
                     lda);
}

inline cublasStatus_t getrf_batched(cublasHandle_t& handle,
                                    int n,
                                    float* A[],
                                    int lda,
                                    int* PivotArray,
                                    int* infoArray,
                                    int batchSize)
{
  return cublasSgetrfBatched(handle, n, A, lda, PivotArray, infoArray, batchSize);
}

inline cublasStatus_t getrf_batched(cublasHandle_t& handle,
                                    int n,
                                    double* A[],
                                    int lda,
                                    int* PivotArray,
                                    int* infoArray,
                                    int batchSize)
{
  return cublasDgetrfBatched(handle, n, A, lda, PivotArray, infoArray, batchSize);
}

inline cublasStatus_t getrf_batched(cublasHandle_t& handle,
                                    int n,
                                    std::complex<float>* A[],
                                    int lda,
                                    int* PivotArray,
                                    int* infoArray,
                                    int batchSize)
{
  return cublasCgetrfBatched(handle, n, castCUDAType(A), lda, PivotArray, infoArray, batchSize);
}

inline cublasStatus_t getrf_batched(cublasHandle_t& handle,
                                    int n,
                                    std::complex<double>* A[],
                                    int lda,
                                    int* PivotArray,
                                    int* infoArray,
                                    int batchSize)
{
  return cublasZgetrfBatched(handle, n, castCUDAType(A), lda, PivotArray, infoArray, batchSize);
}

inline cublasStatus_t getri_batched(cublasHandle_t& handle,
                                    int n,
                                    float* A[],
                                    int lda,
                                    int* PivotArray,
                                    float* C[],
                                    int ldc,
                                    int* infoArray,
                                    int batchSize)
{
  return cublasSgetriBatched(handle, n, A, lda, PivotArray, C, ldc, infoArray, batchSize);
}

inline cublasStatus_t getri_batched(cublasHandle_t& handle,
                                    int n,
                                    double* A[],
                                    int lda,
                                    int* PivotArray,
                                    double* C[],
                                    int ldc,
                                    int* infoArray,
                                    int batchSize)
{
  return cublasDgetriBatched(handle, n, A, lda, PivotArray, C, ldc, infoArray, batchSize);
}

inline cublasStatus_t getri_batched(cublasHandle_t& handle,
                                    int n,
                                    std::complex<float>* A[],
                                    int lda,
                                    int* PivotArray,
                                    std::complex<float>* C[],
                                    int ldc,
                                    int* infoArray,
                                    int batchSize)
{
  return cublasCgetriBatched(handle, n, castCUDAType(A), lda, PivotArray, castCUDAType(C), ldc, infoArray, batchSize);
}

inline cublasStatus_t getri_batched(cublasHandle_t& handle,
                                    int n,
                                    std::complex<double>* A[],
                                    int lda,
                                    int* PivotArray,
                                    std::complex<double>* C[],
                                    int ldc,
                                    int* infoArray,
                                    int batchSize)
{
  return cublasZgetriBatched(handle, n, castCUDAType(A), lda, PivotArray, castCUDAType(C), ldc, infoArray, batchSize);
}

}; // namespace cuBLAS

} // namespace qmcplusplus
#endif // QMCPLUSPLUS_CUBLAS_H
