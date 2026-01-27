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

#include <complex>
#include <iostream>
#include <string>
#include <stdexcept>
#include "config.h"
#ifndef QMC_CUDA2HIP
#include <cublas_v2.h>
#define castNativeType castCUDAType
#else
#include <hipblas/hipblas.h>
#include "Platforms/ROCm/cuda2hip.h"
#include "Platforms/ROCm/hipBLAS.hpp"
#include "Platforms/ROCm/hipblasTypeMapping.hpp"
#define castNativeType casthipblasType
#endif
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
#ifndef QMC_CUDA2HIP
    case CUBLAS_STATUS_LICENSE_ERROR:
      cublas_error = "CUBLAS_STATUS_LICENSE_ERROR";
      break;
#endif
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

inline cublasOperation_t convertOperation(const char trans)
{
  if (trans == 'N' || trans == 'n')
    return CUBLAS_OP_N;
  else if (trans == 'T' || trans == 't')
    return CUBLAS_OP_T;
  else if (trans == 'C' || trans == 'c')
    return CUBLAS_OP_C;
  else
    throw std::runtime_error(
        "cuBLAS::convertOperation trans can only be 'N', 'T', 'C', 'n', 't', 'c'. Input value is " +
        std::string(1, trans));
}

inline cublasStatus_t geam(cublasHandle_t& handle,
                           cublasOperation_t& transa,
                           cublasOperation_t& transb,
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
  const cuDoubleComplex alpha_cu = make_cuDoubleComplex(alpha->real(), alpha->imag());
  const cuDoubleComplex beta_cu  = make_cuDoubleComplex(beta->real(), beta->imag());

  return cublasZgeam(handle, transa, transb, m, n, &alpha_cu, castNativeType(A), lda, &beta_cu,
                     castNativeType(B), ldb, castNativeType(C), ldc);
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
  const cuComplex alpha_cu = make_cuComplex(alpha->real(), alpha->imag());
  const cuComplex beta_cu  = make_cuComplex(beta->real(), beta->imag());

  return cublasCgeam(handle, transa, transb, m, n, &alpha_cu, castNativeType(A), lda, &beta_cu,
                     castNativeType(B), ldb, castNativeType(C), ldc);
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
#undef castNativeType
#endif // QMCPLUSPLUS_CUBLAS_H
