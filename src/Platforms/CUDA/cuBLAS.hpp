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
inline cublasStatus_t gemm(cublasHandle_t& handle,
                           const cublasOperation_t& transa,
                           const cublasOperation_t& transb,
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
  return cublasCgemm(handle, transa, transb, m, n, k, (const cuComplex*)alpha, (const cuComplex*)A, lda,
                     (const cuComplex*)B, ldb, (const cuComplex*)beta, (cuComplex*)C, ldc);
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
  return cublasZgemm(handle, transa, transb, m, n, k, (const cuDoubleComplex*)alpha, (const cuDoubleComplex*)A, lda,
                     (const cuDoubleComplex*)B, ldb, (const cuDoubleComplex*)beta, (cuDoubleComplex*)C, ldc);
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
  return cublasCgemmBatched(handle, transa, transb, m, n, k, (const cuComplex*)alpha, (const cuComplex**)A, lda,
                            (const cuComplex**)B, ldb, (const cuComplex*)beta, (cuComplex**)C, ldc, batchCount);
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
  return cublasZgemmBatched(handle, transa, transb, m, n, k, (const cuDoubleComplex*)alpha, (const cuDoubleComplex**)A,
                            lda, (const cuDoubleComplex**)B, ldb, (const cuDoubleComplex*)beta, (cuDoubleComplex**)C,
                            ldc, batchCount);
}

inline cublasStatus_t matinv_batched(cublasHandle_t& handle,
                                     int n,
                                     const double* A[],
                                     int lda,
                                     double* Ainv[],
                                     int lda_inv,
                                     int* info,
                                     int batchSize)
{
  return cublasDmatinvBatched(handle, n, A, lda, Ainv, lda_inv, info, batchSize);
}

inline cublasStatus_t matinv_batched(cublasHandle_t& handle,
                                     int n,
                                     const float* A[],
                                     int lda,
                                     float* Ainv[],
                                     int lda_inv,
                                     int* info,
                                     int batchSize)
{
  return cublasSmatinvBatched(handle, n, A, lda, Ainv, lda_inv, info, batchSize);
}

inline cublasStatus_t matinv_batched(cublasHandle_t& handle,
                                     int n,
                                     const std::complex<float>* A[],
                                     int lda,
                                     std::complex<float>* Ainv[],
                                     int lda_inv,
                                     int* info,
                                     int batchSize)
{
  return cublasCmatinvBatched(handle, n, reinterpret_cast<const cuFloatComplex**>(A), lda,
                              reinterpret_cast<cuFloatComplex**>(Ainv), lda_inv, info, batchSize);
}

inline cublasStatus_t matinv_batched(cublasHandle_t& handle,
                                     int n,
                                     const std::complex<double>* A[],
                                     int lda,
                                     std::complex<double>* Ainv[],
                                     int lda_inv,
                                     int* info,
                                     int batchSize)
{
  return cublasZmatinvBatched(handle, n, reinterpret_cast<const cuDoubleComplex**>(A), lda,
                              reinterpret_cast<cuDoubleComplex**>(Ainv), lda_inv, info, batchSize);
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
  return cublasCgetrfBatched(handle, n, reinterpret_cast<cuComplex* const *>(A), lda, PivotArray, infoArray, batchSize);
}

inline cublasStatus_t getrf_batched(cublasHandle_t& handle,
                                    int n,
                                    std::complex<double>* A[],
                                    int lda,
                                    int* PivotArray,
                                    int* infoArray,
                                    int batchSize)
{
  return cublasZgetrfBatched(handle, n, reinterpret_cast<cuDoubleComplex* const *>(A), lda, PivotArray, infoArray, batchSize);
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
  return cublasCgetriBatched(handle, n, reinterpret_cast<cuComplex* const *>(A), lda, PivotArray, reinterpret_cast<cuComplex* const *>(C), ldc, infoArray, batchSize);
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
  return cublasZgetriBatched(handle, n, reinterpret_cast<cuDoubleComplex* const *>(A), lda, PivotArray, reinterpret_cast<cuDoubleComplex* const *>(C), ldc, infoArray, batchSize);
}

}; // namespace cuBLAS

} // namespace qmcplusplus
#endif // QMCPLUSPLUS_CUBLAS_H
