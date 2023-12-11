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

#ifndef CUBLAS_FUNCTIONDEFS_H
#define CUBLAS_FUNCTIONDEFS_H

#include <cassert>
#include <cuda_runtime.h>
#include "cublas_v2.h"
#include "AFQMC/Memory/CUDA/cuda_utilities.h"

namespace cublas
{
using qmc_cuda::cublasOperation;

// Level-1
inline cublasStatus_t cublas_copy(cublasHandle_t handle, int n, float* x, int incx, float* y, int incy)
{
  cublasStatus_t success = cublasScopy(handle, n, x, incx, y, incy);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_copy(cublasHandle_t handle, int n, double* x, int incx, double* y, int incy)
{
  cublasStatus_t success = cublasDcopy(handle, n, x, incx, y, incy);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_copy(cublasHandle_t handle,
                                  int n,
                                  std::complex<float>* x,
                                  int incx,
                                  std::complex<float>* y,
                                  int incy)
{
  cublasStatus_t success =
      cublasCcopy(handle, n, reinterpret_cast<cuComplex*>(x), incx, reinterpret_cast<cuComplex*>(y), incy);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_copy(cublasHandle_t handle,
                                  int n,
                                  std::complex<double>* x,
                                  int incx,
                                  std::complex<double>* y,
                                  int incy)
{
  cublasStatus_t success =
      cublasZcopy(handle, n, reinterpret_cast<cuDoubleComplex*>(x), incx, reinterpret_cast<cuDoubleComplex*>(y), incy);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_scal(cublasHandle_t handle, int n, const float alpha, float* x, int incx)
{
  cublasStatus_t success = cublasSscal(handle, n, &alpha, x, incx);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_scal(cublasHandle_t handle, int n, const double alpha, double* x, int incx)
{
  cublasStatus_t success = cublasDscal(handle, n, &alpha, x, incx);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_scal(cublasHandle_t handle,
                                  int n,
                                  const std::complex<float> alpha,
                                  std::complex<float>* x,
                                  int incx)
{
  cublasStatus_t success =
      cublasCscal(handle, n, reinterpret_cast<cuComplex const*>(&alpha), reinterpret_cast<cuComplex*>(x), incx);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_scal(cublasHandle_t handle,
                                  int n,
                                  const std::complex<double> alpha,
                                  std::complex<double>* x,
                                  int incx)
{
  cublasStatus_t success = cublasZscal(handle, n, reinterpret_cast<cuDoubleComplex const*>(&alpha),
                                      reinterpret_cast<cuDoubleComplex*>(x), incx);
  cudaDeviceSynchronize();
  return success;
}

inline float cublas_dot(cublasHandle_t handle, int n, const float* x, int incx, const float* y, int incy)
{
  float result;
  cublasStatus_t success = cublasSdot(handle, n, x, incx, y, incy, &result);
  cudaDeviceSynchronize();
  if (CUBLAS_STATUS_SUCCESS != success)
    throw std::runtime_error("Error: cublas_dot returned error code.");
  return result;
}

inline double cublas_dot(cublasHandle_t handle, int n, const double* x, int incx, const double* y, int incy)
{
  double result;
  cublasStatus_t success = cublasDdot(handle, n, x, incx, y, incy, &result);
  cudaDeviceSynchronize();
  if (CUBLAS_STATUS_SUCCESS != success)
    throw std::runtime_error("Error: cublas_dot returned error code.");
  return result;
}

inline std::complex<float> cublas_dot(cublasHandle_t handle,
                                      int n,
                                      const std::complex<float>* x,
                                      int incx,
                                      const std::complex<float>* y,
                                      int incy)
{
  std::complex<float> result;
  cublasStatus_t success =
      cublasCdotu(handle, n, reinterpret_cast<cuComplex const*>(x), incx, reinterpret_cast<cuComplex const*>(y), incy,
                  reinterpret_cast<cuComplex*>(&result));
  cudaDeviceSynchronize();
  if (CUBLAS_STATUS_SUCCESS != success)
    throw std::runtime_error("Error: cublas_dot returned error code.");
  return result;
}

inline std::complex<double> cublas_dot(cublasHandle_t handle,
                                       int n,
                                       const std::complex<double>* x,
                                       int incx,
                                       const std::complex<double>* y,
                                       int incy)
{
  std::complex<double> result;
  cublasStatus_t success =
      cublasZdotu(handle, n, reinterpret_cast<cuDoubleComplex const*>(x), incx,
                  reinterpret_cast<cuDoubleComplex const*>(y), incy, reinterpret_cast<cuDoubleComplex*>(&result));
  cudaDeviceSynchronize();
  if (CUBLAS_STATUS_SUCCESS != success)
    throw std::runtime_error("Error: cublas_dot returned error code.");
  return result;
}

inline std::complex<double> cublas_dot(cublasHandle_t handle,
                                       int n,
                                       const double* x,
                                       int incx,
                                       const std::complex<double>* y,
                                       int incy)
{
  int incy_         = 2 * incy;
  const double* y_  = reinterpret_cast<const double*>(y);
  const double* y1_ = y_ + 1;
  double resR, resI;
  cublasStatus_t success = cublasDdot(handle, n, x, incx, y_, incy_, &resR);
  cudaDeviceSynchronize();
  if (CUBLAS_STATUS_SUCCESS != success)
    throw std::runtime_error("Error: cublas_dot returned error code.");
  success = cublasDdot(handle, n, x, incx, y1_, incy_, &resI);
  cudaDeviceSynchronize();
  if (CUBLAS_STATUS_SUCCESS != success)
    throw std::runtime_error("Error: cublas_dot returned error code.");
  return std::complex<double>{resR, resI};
}

inline std::complex<double> cublas_dot(cublasHandle_t handle,
                                       int n,
                                       const std::complex<double>* x,
                                       int incx,
                                       const double* y,
                                       int incy)
{
  int incx_         = 2 * incx;
  const double* x_  = reinterpret_cast<const double*>(x);
  const double* x1_ = x_ + 1;
  double resR, resI;
  cublasStatus_t success = cublasDdot(handle, n, x_, incx_, y, incy, &resR);
  cudaDeviceSynchronize();
  if (CUBLAS_STATUS_SUCCESS != success)
    throw std::runtime_error("Error: cublas_dot returned error code.");
  success = cublasDdot(handle, n, x1_, incx_, y, incy, &resI);
  cudaDeviceSynchronize();
  if (CUBLAS_STATUS_SUCCESS != success)
    throw std::runtime_error("Error: cublas_dot returned error code.");
  return std::complex<double>{resR, resI};
}

inline cublasStatus_t cublas_axpy(cublasHandle_t handle,
                                  int n,
                                  const float alpha,
                                  const float* x,
                                  int incx,
                                  float* y,
                                  int incy)
{
  cublasStatus_t success = cublasSaxpy(handle, n, &alpha, x, incx, y, incy);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_axpy(cublasHandle_t handle,
                                  int n,
                                  const double alpha,
                                  const double* x,
                                  int incx,
                                  double* y,
                                  int incy)
{
  cublasStatus_t success = cublasDaxpy(handle, n, &alpha, x, incx, y, incy);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_axpy(cublasHandle_t handle,
                                  int n,
                                  const std::complex<float> alpha,
                                  const std::complex<float>* x,
                                  int incx,
                                  std::complex<float>* y,
                                  int incy)
{
  cublasStatus_t success =
      cublasCaxpy(handle, n, reinterpret_cast<cuComplex const*>(&alpha), reinterpret_cast<cuComplex const*>(x), incx,
                  reinterpret_cast<cuComplex*>(y), incy);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_axpy(cublasHandle_t handle,
                                  int n,
                                  const std::complex<double> alpha,
                                  const std::complex<double>* x,
                                  int incx,
                                  std::complex<double>* y,
                                  int incy)
{
  cublasStatus_t success =
      cublasZaxpy(handle, n, reinterpret_cast<cuDoubleComplex const*>(&alpha),
                  reinterpret_cast<cuDoubleComplex const*>(x), incx, reinterpret_cast<cuDoubleComplex*>(y), incy);
  cudaDeviceSynchronize();
  return success;
}

// Level-2
inline cublasStatus_t cublas_gemv(cublasHandle_t handle,
                                  char Atrans,
                                  int M,
                                  int N,
                                  const float alpha,
                                  const float* A,
                                  int lda,
                                  const float* x,
                                  int incx,
                                  const float beta,
                                  float* y,
                                  int incy)
{
  cublasStatus_t success = cublasSgemv(handle, cublasOperation(Atrans), M, N, &alpha, A, lda, x, incx, &beta, y, incy);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_gemv(cublasHandle_t handle,
                                  char Atrans,
                                  int M,
                                  int N,
                                  const double alpha,
                                  const double* A,
                                  int lda,
                                  const double* x,
                                  int incx,
                                  const double beta,
                                  double* y,
                                  int incy)
{
  cublasStatus_t success = cublasDgemv(handle, cublasOperation(Atrans), M, N, &alpha, A, lda, x, incx, &beta, y, incy);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_gemv(cublasHandle_t handle,
                                  char Atrans,
                                  int M,
                                  int N,
                                  const std::complex<float> alpha,
                                  const std::complex<float>* A,
                                  int lda,
                                  const std::complex<float>* x,
                                  int incx,
                                  const std::complex<float> beta,
                                  std::complex<float>* y,
                                  int incy)
{
  cublasStatus_t success =
      cublasCgemv(handle, cublasOperation(Atrans), M, N, reinterpret_cast<cuComplex const*>(&alpha),
                  reinterpret_cast<cuComplex const*>(A), lda, reinterpret_cast<cuComplex const*>(x), incx,
                  reinterpret_cast<cuComplex const*>(&beta), reinterpret_cast<cuComplex*>(y), incy);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_gemv(cublasHandle_t handle,
                                  char Atrans,
                                  int M,
                                  int N,
                                  const std::complex<double> alpha,
                                  const std::complex<double>* A,
                                  int lda,
                                  const std::complex<double>* x,
                                  int incx,
                                  const std::complex<double> beta,
                                  std::complex<double>* y,
                                  int incy)
{
  cublasStatus_t success =
      cublasZgemv(handle, cublasOperation(Atrans), M, N, reinterpret_cast<cuDoubleComplex const*>(&alpha),
                  reinterpret_cast<cuDoubleComplex const*>(A), lda, reinterpret_cast<cuDoubleComplex const*>(x), incx,
                  reinterpret_cast<cuDoubleComplex const*>(&beta), reinterpret_cast<cuDoubleComplex*>(y), incy);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_gemv(cublasHandle_t handle,
                                  char Atrans,
                                  int M,
                                  int N,
                                  const float alpha,
                                  const float* A,
                                  int lda,
                                  const std::complex<float>* x,
                                  int incx,
                                  const float beta,
                                  std::complex<float>* y,
                                  int incy)
{
  cublasStatus_t success = CUBLAS_STATUS_SUCCESS;
  char Nt('N');
  char Tt('T');
  if (Atrans == 'n' || Atrans == 'N')
    success =
        cublasSgemm(handle, cublasOperation(Nt), cublasOperation(Tt), 2, M, N, &alpha,
                    reinterpret_cast<float const*>(x), 2 * incx, A, lda, &beta, reinterpret_cast<float*>(y), 2 * incy);
  else if (Atrans == 't' || Atrans == 'T')
    success =
        cublasSgemm(handle, cublasOperation(Nt), cublasOperation(Nt), 2, N, M, &alpha,
                    reinterpret_cast<float const*>(x), 2 * incx, A, lda, &beta, reinterpret_cast<float*>(y), 2 * incy);
  else
    assert(0);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_gemv(cublasHandle_t handle,
                                  char Atrans,
                                  int M,
                                  int N,
                                  const double alpha,
                                  const double* A,
                                  int lda,
                                  const std::complex<double>* x,
                                  int incx,
                                  const double beta,
                                  std::complex<double>* y,
                                  int incy)
{
  cublasStatus_t success = CUBLAS_STATUS_SUCCESS;
  char Nt('N');
  char Tt('T');
  if (Atrans == 'n' || Atrans == 'N')
    success = cublasDgemm(handle, cublasOperation(Nt), cublasOperation(Tt), 2, M, N, &alpha,
                         reinterpret_cast<double const*>(x), 2 * incx, A, lda, &beta, reinterpret_cast<double*>(y),
                         2 * incy);
  else if (Atrans == 't' || Atrans == 'T')
    success = cublasDgemm(handle, cublasOperation(Nt), cublasOperation(Nt), 2, N, M, &alpha,
                         reinterpret_cast<double const*>(x), 2 * incx, A, lda, &beta, reinterpret_cast<double*>(y),
                         2 * incy);
  else
    assert(0);
  cudaDeviceSynchronize();
  return success;
}


// Level-3
inline cublasStatus_t cublas_gemm(cublasHandle_t handle,
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
  cublasStatus_t success = cublasSgemm(handle, cublasOperation(Atrans), cublasOperation(Btrans), M, N, K, &alpha, A, lda,
                                      B, ldb, &beta, C, ldc);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_gemm(cublasHandle_t handle,
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
  cublasStatus_t success = cublasDgemm(handle, cublasOperation(Atrans), cublasOperation(Btrans), M, N, K, &alpha, A, lda,
                                      B, ldb, &beta, C, ldc);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_gemm(cublasHandle_t handle,
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
  cublasStatus_t success = cublasCgemm(handle, cublasOperation(Atrans), cublasOperation(Btrans), M, N, K,
                                      reinterpret_cast<cuComplex const*>(&alpha), reinterpret_cast<cuComplex const*>(A),
                                      lda, reinterpret_cast<cuComplex const*>(B), ldb,
                                      reinterpret_cast<cuComplex const*>(&beta), reinterpret_cast<cuComplex*>(C), ldc);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_gemm(cublasHandle_t handle,
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
      cublasZgemm(handle, cublasOperation(Atrans), cublasOperation(Btrans), M, N, K,
                  reinterpret_cast<cuDoubleComplex const*>(&alpha), reinterpret_cast<cuDoubleComplex const*>(A), lda,
                  reinterpret_cast<cuDoubleComplex const*>(B), ldb, reinterpret_cast<cuDoubleComplex const*>(&beta),
                  reinterpret_cast<cuDoubleComplex*>(C), ldc);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_gemm(cublasHandle_t handle,
                                  char Atrans,
                                  char Btrans,
                                  int M,
                                  int N,
                                  int K,
                                  const float alpha,
                                  const std::complex<float>* A,
                                  int lda,
                                  const float* B,
                                  int ldb,
                                  const float beta,
                                  std::complex<float>* C,
                                  int ldc)
{
  assert(Atrans == 'n' || Atrans == 'N');
  cublasStatus_t success =
      cublasSgemm(handle, cublasOperation(Atrans), cublasOperation(Btrans), 2 * M, N, K, &alpha,
                  reinterpret_cast<float const*>(A), 2 * lda, B, ldb, &beta, reinterpret_cast<float*>(C), 2 * ldc);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_gemm(cublasHandle_t handle,
                                  char Atrans,
                                  char Btrans,
                                  int M,
                                  int N,
                                  int K,
                                  double alpha,
                                  const std::complex<double>* A,
                                  int lda,
                                  const double* B,
                                  int ldb,
                                  double beta,
                                  std::complex<double>* C,
                                  int ldc)
{
  assert(Atrans == 'n' || Atrans == 'N');
  cublasStatus_t success =
      cublasDgemm(handle, cublasOperation(Atrans), cublasOperation(Btrans), 2 * M, N, K, &alpha,
                  reinterpret_cast<double const*>(A), 2 * lda, B, ldb, &beta, reinterpret_cast<double*>(C), 2 * ldc);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_gemm(cublasHandle_t handle,
                                  char Atrans,
                                  char Btrans,
                                  int M,
                                  int N,
                                  int K,
                                  const cuDoubleComplex alpha,
                                  const cuDoubleComplex* A,
                                  int lda,
                                  const cuDoubleComplex* B,
                                  int ldb,
                                  const cuDoubleComplex beta,
                                  cuDoubleComplex* C,
                                  int ldc)
{
  cublasStatus_t success = cublasZgemm(handle, cublasOperation(Atrans), cublasOperation(Btrans), M, N, K, &alpha, A, lda,
                                      B, ldb, &beta, C, ldc);
  cudaDeviceSynchronize();
  return success;
}

// Extensions
inline cublasStatus_t cublas_getrfBatched(cublasHandle_t handle,
                                          int n,
                                          float** Aarray,
                                          int lda,
                                          int* PivotArray,
                                          int* infoArray,
                                          int batchSize)
{
  cublasStatus_t success = cublasSgetrfBatched(handle, n, Aarray, lda, PivotArray, infoArray, batchSize);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_getrfBatched(cublasHandle_t handle,
                                          int n,
                                          double** Aarray,
                                          int lda,
                                          int* PivotArray,
                                          int ldc,
                                          int* infoArray,
                                          int batchSize)
{
  cublasStatus_t success = cublasDgetrfBatched(handle, n, Aarray, lda, PivotArray, infoArray, batchSize);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_getrfBatched(cublasHandle_t handle,
                                          int n,
                                          std::complex<double>** Aarray,
                                          int lda,
                                          int* PivotArray,
                                          int* infoArray,
                                          int batchSize)
{
  cublasStatus_t success = cublasZgetrfBatched(handle, n, reinterpret_cast<cuDoubleComplex* const*>(Aarray), lda,
                                              PivotArray, infoArray, batchSize);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_getrfBatched(cublasHandle_t handle,
                                          int n,
                                          std::complex<float>** Aarray,
                                          int lda,
                                          int* PivotArray,
                                          int* infoArray,
                                          int batchSize)
{
  cublasStatus_t success = cublasCgetrfBatched(handle, n, reinterpret_cast<cuComplex* const*>(Aarray), lda, PivotArray,
                                              infoArray, batchSize);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_getriBatched(cublasHandle_t handle,
                                          int n,
                                          float** Aarray,
                                          int lda,
                                          int* PivotArray,
                                          float** Carray,
                                          int ldc,
                                          int* infoArray,
                                          int batchSize)
{
  cublasStatus_t success = cublasSgetriBatched(handle, n, Aarray, lda, PivotArray, Carray, ldc, infoArray, batchSize);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_getriBatched(cublasHandle_t handle,
                                          int n,
                                          double** Aarray,
                                          int lda,
                                          int* PivotArray,
                                          double** Carray,
                                          int ldc,
                                          int* infoArray,
                                          int batchSize)
{
  cublasStatus_t success = cublasDgetriBatched(handle, n, Aarray, lda, PivotArray, Carray, ldc, infoArray, batchSize);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_getriBatched(cublasHandle_t handle,
                                          int n,
                                          std::complex<double>** Aarray,
                                          int lda,
                                          int* PivotArray,
                                          std::complex<double>** Carray,
                                          int ldc,
                                          int* infoArray,
                                          int batchSize)
{
  cublasStatus_t success =
      cublasZgetriBatched(handle, n, reinterpret_cast<const cuDoubleComplex* const*>(Aarray), lda, PivotArray,
                          reinterpret_cast<cuDoubleComplex* const*>(Carray), ldc, infoArray, batchSize);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_getriBatched(cublasHandle_t handle,
                                          int n,
                                          std::complex<float>** Aarray,
                                          int lda,
                                          int* PivotArray,
                                          std::complex<float>** Carray,
                                          int ldc,
                                          int* infoArray,
                                          int batchSize)
{
  cublasStatus_t success =
      cublasCgetriBatched(handle, n, reinterpret_cast<const cuComplex* const*>(Aarray), lda, PivotArray,
                          reinterpret_cast<cuComplex* const*>(Carray), ldc, infoArray, batchSize);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_matinvBatched(cublasHandle_t handle,
                                           int n,
                                           float** Aarray,
                                           int lda,
                                           float** Carray,
                                           int ldc,
                                           int* infoArray,
                                           int batchSize)
{
  cublasStatus_t success = cublasSmatinvBatched(handle, n, Aarray, lda, Carray, ldc, infoArray, batchSize);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_matinvBatched(cublasHandle_t handle,
                                           int n,
                                           double** Aarray,
                                           int lda,
                                           double** Carray,
                                           int ldc,
                                           int* infoArray,
                                           int batchSize)
{
  cublasStatus_t success = cublasDmatinvBatched(handle, n, Aarray, lda, Carray, ldc, infoArray, batchSize);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_matinvBatched(cublasHandle_t handle,
                                           int n,
                                           std::complex<float>** Aarray,
                                           int lda,
                                           std::complex<float>** Carray,
                                           int ldc,
                                           int* infoArray,
                                           int batchSize)
{
  cublasStatus_t success = cublasCmatinvBatched(handle, n, reinterpret_cast<const cuComplex* const*>(Aarray), lda,
                                               reinterpret_cast<cuComplex**>(Carray), ldc, infoArray, batchSize);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_matinvBatched(cublasHandle_t handle,
                                           int n,
                                           std::complex<double>** Aarray,
                                           int lda,
                                           std::complex<double>** Carray,
                                           int ldc,
                                           int* infoArray,
                                           int batchSize)
{
  cublasStatus_t success = cublasZmatinvBatched(handle, n, reinterpret_cast<const cuDoubleComplex* const*>(Aarray), lda,
                                               reinterpret_cast<cuDoubleComplex**>(Carray), ldc, infoArray, batchSize);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_geam(cublasHandle_t handle,
                                  char Atrans,
                                  char Btrans,
                                  int M,
                                  int N,
                                  float const alpha,
                                  float const* A,
                                  int lda,
                                  float const beta,
                                  float const* B,
                                  int ldb,
                                  float* C,
                                  int ldc)
{
  cublasStatus_t success = cublasSgeam(handle, cublasOperation(Atrans), cublasOperation(Btrans), M, N, &alpha, A, lda,
                                      &beta, B, ldb, C, ldc);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_geam(cublasHandle_t handle,
                                  char Atrans,
                                  char Btrans,
                                  int M,
                                  int N,
                                  double const alpha,
                                  double const* A,
                                  int lda,
                                  double const beta,
                                  double const* B,
                                  int ldb,
                                  double* C,
                                  int ldc)
{
  cublasStatus_t success = cublasDgeam(handle, cublasOperation(Atrans), cublasOperation(Btrans), M, N, &alpha, A, lda,
                                      &beta, B, ldb, C, ldc);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_geam(cublasHandle_t handle,
                                  char Atrans,
                                  char Btrans,
                                  int M,
                                  int N,
                                  std::complex<float> const alpha,
                                  std::complex<float> const* A,
                                  int lda,
                                  std::complex<float> const beta,
                                  std::complex<float> const* B,
                                  int ldb,
                                  std::complex<float>* C,
                                  int ldc)
{
  cublasStatus_t success = cublasCgeam(handle, cublasOperation(Atrans), cublasOperation(Btrans), M, N,
                                      reinterpret_cast<cuComplex const*>(&alpha), reinterpret_cast<cuComplex const*>(A),
                                      lda, reinterpret_cast<cuComplex const*>(&beta),
                                      reinterpret_cast<cuComplex const*>(B), ldb, reinterpret_cast<cuComplex*>(C), ldc);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_geam(cublasHandle_t handle,
                                  char Atrans,
                                  char Btrans,
                                  int M,
                                  int N,
                                  std::complex<double> const alpha,
                                  std::complex<double> const* A,
                                  int lda,
                                  std::complex<double> const beta,
                                  std::complex<double> const* B,
                                  int ldb,
                                  std::complex<double>* C,
                                  int ldc)
{
  cublasStatus_t success =
      cublasZgeam(handle, cublasOperation(Atrans), cublasOperation(Btrans), M, N,
                  reinterpret_cast<cuDoubleComplex const*>(&alpha), reinterpret_cast<cuDoubleComplex const*>(A), lda,
                  reinterpret_cast<cuDoubleComplex const*>(&beta), reinterpret_cast<cuDoubleComplex const*>(B), ldb,
                  reinterpret_cast<cuDoubleComplex*>(C), ldc);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_gemmStridedBatched(cublasHandle_t handle,
                                                char Atrans,
                                                char Btrans,
                                                int M,
                                                int N,
                                                int K,
                                                const float alpha,
                                                const float* A,
                                                int lda,
                                                int strideA,
                                                const float* B,
                                                int ldb,
                                                int strideB,
                                                const float beta,
                                                float* C,
                                                int ldc,
                                                int strideC,
                                                int batchSize)
{
  cublasStatus_t success =
      cublasSgemmStridedBatched(handle, cublasOperation(Atrans), cublasOperation(Btrans), M, N, K, &alpha, A, lda,
                                strideA, B, ldb, strideB, &beta, C, ldc, strideC, batchSize);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_gemmStridedBatched(cublasHandle_t handle,
                                                char Atrans,
                                                char Btrans,
                                                int M,
                                                int N,
                                                int K,
                                                const double alpha,
                                                const double* A,
                                                int lda,
                                                int strideA,
                                                const double* B,
                                                int ldb,
                                                int strideB,
                                                const double beta,
                                                double* C,
                                                int ldc,
                                                int strideC,
                                                int batchSize)
{
  cublasStatus_t success =
      cublasDgemmStridedBatched(handle, cublasOperation(Atrans), cublasOperation(Btrans), M, N, K, &alpha, A, lda,
                                strideA, B, ldb, strideB, &beta, C, ldc, strideC, batchSize);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_gemmStridedBatched(cublasHandle_t handle,
                                                char Atrans,
                                                char Btrans,
                                                int M,
                                                int N,
                                                int K,
                                                const std::complex<float> alpha,
                                                const std::complex<float>* A,
                                                int lda,
                                                int strideA,
                                                const std::complex<float>* B,
                                                int ldb,
                                                int strideB,
                                                const std::complex<float> beta,
                                                std::complex<float>* C,
                                                int ldc,
                                                int strideC,
                                                int batchSize)
{
  cublasStatus_t success =
      cublasCgemmStridedBatched(handle, cublasOperation(Atrans), cublasOperation(Btrans), M, N, K,
                                reinterpret_cast<cuComplex const*>(&alpha), reinterpret_cast<cuComplex const*>(A), lda,
                                strideA, reinterpret_cast<cuComplex const*>(B), ldb, strideB,
                                reinterpret_cast<cuComplex const*>(&beta), reinterpret_cast<cuComplex*>(C), ldc,
                                strideC, batchSize);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_gemmStridedBatched(cublasHandle_t handle,
                                                char Atrans,
                                                char Btrans,
                                                int M,
                                                int N,
                                                int K,
                                                const std::complex<double> alpha,
                                                const std::complex<double>* A,
                                                int lda,
                                                int strideA,
                                                const std::complex<double>* B,
                                                int ldb,
                                                int strideB,
                                                const std::complex<double> beta,
                                                std::complex<double>* C,
                                                int ldc,
                                                int strideC,
                                                int batchSize)
{
  cublasStatus_t success = cublasZgemmStridedBatched(handle, cublasOperation(Atrans), cublasOperation(Btrans), M, N, K,
                                                    reinterpret_cast<cuDoubleComplex const*>(&alpha),
                                                    reinterpret_cast<cuDoubleComplex const*>(A), lda, strideA,
                                                    reinterpret_cast<cuDoubleComplex const*>(B), ldb, strideB,
                                                    reinterpret_cast<cuDoubleComplex const*>(&beta),
                                                    reinterpret_cast<cuDoubleComplex*>(C), ldc, strideC, batchSize);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_gemmBatched(cublasHandle_t handle,
                                         char Atrans,
                                         char Btrans,
                                         int M,
                                         int N,
                                         int K,
                                         float alpha,
                                         float** A,
                                         int lda,
                                         float** B,
                                         int ldb,
                                         float beta,
                                         float** C,
                                         int ldc,
                                         int batchSize)
{
  cublasStatus_t success = cublasSgemmBatched(handle, cublasOperation(Atrans), cublasOperation(Btrans), M, N, K, &alpha,
                                             A, lda, B, ldb, &beta, C, ldc, batchSize);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_gemmBatched(cublasHandle_t handle,
                                         char Atrans,
                                         char Btrans,
                                         int M,
                                         int N,
                                         int K,
                                         double alpha,
                                         double** A,
                                         int lda,
                                         double** B,
                                         int ldb,
                                         double beta,
                                         double** C,
                                         int ldc,
                                         int batchSize)
{
  cublasStatus_t success = cublasDgemmBatched(handle, cublasOperation(Atrans), cublasOperation(Btrans), M, N, K, &alpha,
                                             A, lda, B, ldb, &beta, C, ldc, batchSize);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_gemmBatched(cublasHandle_t handle,
                                         char Atrans,
                                         char Btrans,
                                         int M,
                                         int N,
                                         int K,
                                         std::complex<float> alpha,
                                         std::complex<float>** A,
                                         int lda,
                                         std::complex<float>** B,
                                         int ldb,
                                         std::complex<float> beta,
                                         std::complex<float>** C,
                                         int ldc,
                                         int batchSize)
{
  cublasStatus_t success =
      cublasCgemmBatched(handle, cublasOperation(Atrans), cublasOperation(Btrans), M, N, K,
                         reinterpret_cast<cuComplex*>(&alpha), reinterpret_cast<cuComplex**>(A), lda,
                         reinterpret_cast<cuComplex**>(B), ldb, reinterpret_cast<cuComplex*>(&beta),
                         reinterpret_cast<cuComplex**>(C), ldc, batchSize);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_gemmBatched(cublasHandle_t handle,
                                         char Atrans,
                                         char Btrans,
                                         int M,
                                         int N,
                                         int K,
                                         std::complex<double> alpha,
                                         std::complex<double>** A,
                                         int lda,
                                         std::complex<double>** B,
                                         int ldb,
                                         std::complex<double> beta,
                                         std::complex<double>** C,
                                         int ldc,
                                         int batchSize)
{
  cublasStatus_t success =
      cublasZgemmBatched(handle, cublasOperation(Atrans), cublasOperation(Btrans), M, N, K,
                         reinterpret_cast<cuDoubleComplex*>(&alpha), reinterpret_cast<cuDoubleComplex**>(A), lda,
                         reinterpret_cast<cuDoubleComplex**>(B), ldb, reinterpret_cast<cuDoubleComplex*>(&beta),
                         reinterpret_cast<cuDoubleComplex**>(C), ldc, batchSize);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_gemmBatched(cublasHandle_t handle,
                                         char Atrans,
                                         char Btrans,
                                         int M,
                                         int N,
                                         int K,
                                         float alpha,
                                         std::complex<float>** A,
                                         int lda,
                                         float** B,
                                         int ldb,
                                         float beta,
                                         std::complex<float>** C,
                                         int ldc,
                                         int batchSize)
{
  cublasStatus_t success = cublasSgemmBatched(handle, cublasOperation(Atrans), cublasOperation(Btrans), 2 * M, N, K,
                                             &alpha, reinterpret_cast<float**>(A), 2 * lda, B, ldb, &beta,
                                             reinterpret_cast<float**>(C), 2 * ldc, batchSize);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_gemmBatched(cublasHandle_t handle,
                                         char Atrans,
                                         char Btrans,
                                         int M,
                                         int N,
                                         int K,
                                         double alpha,
                                         std::complex<double>** A,
                                         int lda,
                                         double** B,
                                         int ldb,
                                         double beta,
                                         std::complex<double>** C,
                                         int ldc,
                                         int batchSize)
{
  cublasStatus_t success = cublasDgemmBatched(handle, cublasOperation(Atrans), cublasOperation(Btrans), 2 * M, N, K,
                                             &alpha, reinterpret_cast<double**>(A), 2 * lda, B, ldb, &beta,
                                             reinterpret_cast<double**>(C), 2 * ldc, batchSize);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_geqrfBatched(cublasHandle_t handle,
                                          int m,
                                          int n,
                                          double** Aarray,
                                          int lda,
                                          double** TauArray,
                                          int* info,
                                          int batchSize)
{
  cublasStatus_t success = cublasDgeqrfBatched(handle, m, n, Aarray, lda, TauArray, info, batchSize);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_geqrfBatched(cublasHandle_t handle,
                                          int m,
                                          int n,
                                          float** Aarray,
                                          int lda,
                                          float** TauArray,
                                          int* info,
                                          int batchSize)
{
  cublasStatus_t success = cublasSgeqrfBatched(handle, m, n, Aarray, lda, TauArray, info, batchSize);
  cudaDeviceSynchronize();
  return success;
}


inline cublasStatus_t cublas_geqrfBatched(cublasHandle_t handle,
                                          int m,
                                          int n,
                                          std::complex<double>** Aarray,
                                          int lda,
                                          std::complex<double>** TauArray,
                                          int* info,
                                          int batchSize)
{
  cublasStatus_t success = cublasZgeqrfBatched(handle, m, n, reinterpret_cast<cuDoubleComplex**>(Aarray), lda,
                                              reinterpret_cast<cuDoubleComplex**>(TauArray), info, batchSize);
  cudaDeviceSynchronize();
  return success;
}

inline cublasStatus_t cublas_geqrfBatched(cublasHandle_t handle,
                                          int m,
                                          int n,
                                          std::complex<float>** Aarray,
                                          int lda,
                                          std::complex<float>** TauArray,
                                          int* info,
                                          int batchSize)
{
  cublasStatus_t success = cublasCgeqrfBatched(handle, m, n, reinterpret_cast<cuComplex**>(Aarray), lda,
                                              reinterpret_cast<cuComplex**>(TauArray), info, batchSize);
  cudaDeviceSynchronize();
  return success;
}

} // namespace cublas

#endif
