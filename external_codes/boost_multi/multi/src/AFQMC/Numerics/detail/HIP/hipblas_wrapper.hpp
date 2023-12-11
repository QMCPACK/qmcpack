///////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Fionn Malone, malone14@llnl.gov, LLNL
//
// File created by: Fionn Malone, malone14@llnl.gov, LLNL
////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_HIPBLAS_FUNCTIONDEFS_H
#define AFQMC_HIPBLAS_FUNCTIONDEFS_H

#include <cassert>
#include <hip/hip_runtime.h>
#include <hipblas/hipblas.h>
#include <rocsolver/rocsolver.h>
#include "AFQMC/Memory/HIP/hip_utilities.h"

namespace hipblas
{
using qmc_hip::hipblasOperation;

// Level-1
inline hipblasStatus_t hipblas_copy(hipblasHandle_t handle, int n, float* x, int incx, float* y, int incy)
{
  hipblasStatus_t success = hipblasScopy(handle, n, x, incx, y, incy);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_copy(hipblasHandle_t handle, int n, double* x, int incx, double* y, int incy)
{
  hipblasStatus_t success = hipblasDcopy(handle, n, x, incx, y, incy);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_copy(hipblasHandle_t handle,
                                    int n,
                                    std::complex<float>* x,
                                    int incx,
                                    std::complex<float>* y,
                                    int incy)
{
  hipblasStatus_t success =
      hipblasCcopy(handle, n, reinterpret_cast<hipblasComplex*>(x), incx, reinterpret_cast<hipblasComplex*>(y), incy);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_copy(hipblasHandle_t handle,
                                    int n,
                                    std::complex<double>* x,
                                    int incx,
                                    std::complex<double>* y,
                                    int incy)
{
  hipblasStatus_t success = hipblasZcopy(handle, n, reinterpret_cast<hipblasDoubleComplex*>(x), incx,
                                         reinterpret_cast<hipblasDoubleComplex*>(y), incy);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_scal(hipblasHandle_t handle, int n, const float alpha, float* x, int incx)
{
  hipblasStatus_t success = hipblasSscal(handle, n, &alpha, x, incx);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_scal(hipblasHandle_t handle, int n, const double alpha, double* x, int incx)
{
  hipblasStatus_t success = hipblasDscal(handle, n, &alpha, x, incx);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_scal(hipblasHandle_t handle,
                                    int n,
                                    const std::complex<float> alpha,
                                    std::complex<float>* x,
                                    int incx)
{
  hipblasStatus_t success = hipblasCscal(handle, n, reinterpret_cast<hipblasComplex const*>(&alpha),
                                         reinterpret_cast<hipblasComplex*>(x), incx);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_scal(hipblasHandle_t handle,
                                    int n,
                                    const std::complex<double> alpha,
                                    std::complex<double>* x,
                                    int incx)
{
  hipblasStatus_t success = hipblasZscal(handle, n, reinterpret_cast<hipblasDoubleComplex const*>(&alpha),
                                         reinterpret_cast<hipblasDoubleComplex*>(x), incx);
  hipDeviceSynchronize();
  return success;
}

inline float hipblas_dot(hipblasHandle_t handle, int n, const float* x, int incx, const float* y, int incy)
{
  float result;
  hipblasStatus_t success = hipblasSdot(handle, n, x, incx, y, incy, &result);
  hipDeviceSynchronize();
  if (HIPBLAS_STATUS_SUCCESS != success)
    throw std::runtime_error("Error: hipblas_dot returned error code.");
  return result;
}

inline double hipblas_dot(hipblasHandle_t handle, int n, const double* x, int incx, const double* y, int incy)
{
  double result;
  hipblasStatus_t success = hipblasDdot(handle, n, x, incx, y, incy, &result);
  hipDeviceSynchronize();
  if (HIPBLAS_STATUS_SUCCESS != success)
    throw std::runtime_error("Error: hipblas_dot returned error code.");
  return result;
}

inline std::complex<float> hipblas_dot(hipblasHandle_t handle,
                                       int n,
                                       const std::complex<float>* x,
                                       int incx,
                                       const std::complex<float>* y,
                                       int incy)
{
  std::complex<float> result;
  hipblasStatus_t success =
      hipblasCdotu(handle, n, reinterpret_cast<hipblasComplex const*>(x), incx,
                   reinterpret_cast<hipblasComplex const*>(y), incy, reinterpret_cast<hipblasComplex*>(&result));
  hipDeviceSynchronize();
  if (HIPBLAS_STATUS_SUCCESS != success)
    throw std::runtime_error("Error: hipblas_dot returned error code.");
  return result;
}

inline std::complex<double> hipblas_dot(hipblasHandle_t handle,
                                        int n,
                                        const std::complex<double>* x,
                                        int incx,
                                        const std::complex<double>* y,
                                        int incy)
{
  std::complex<double> result;
  hipblasStatus_t success = hipblasZdotu(handle, n, reinterpret_cast<hipblasDoubleComplex const*>(x), incx,
                                         reinterpret_cast<hipblasDoubleComplex const*>(y), incy,
                                         reinterpret_cast<hipblasDoubleComplex*>(&result));
  hipDeviceSynchronize();
  if (HIPBLAS_STATUS_SUCCESS != success)
    throw std::runtime_error("Error: hipblas_dot returned error code.");
  return result;
}

inline std::complex<double> hipblas_dot(hipblasHandle_t handle,
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
  hipblasStatus_t success = hipblasDdot(handle, n, x, incx, y_, incy_, &resR);
  hipDeviceSynchronize();
  if (HIPBLAS_STATUS_SUCCESS != success)
    throw std::runtime_error("Error: hipblas_dot returned error code.");
  success = hipblasDdot(handle, n, x, incx, y1_, incy_, &resI);
  hipDeviceSynchronize();
  if (HIPBLAS_STATUS_SUCCESS != success)
    throw std::runtime_error("Error: hipblas_dot returned error code.");
  return std::complex<double>{resR, resI};
}

inline std::complex<double> hipblas_dot(hipblasHandle_t handle,
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
  hipblasStatus_t success = hipblasDdot(handle, n, x_, incx_, y, incy, &resR);
  hipDeviceSynchronize();
  if (HIPBLAS_STATUS_SUCCESS != success)
    throw std::runtime_error("Error: hipblas_dot returned error code.");
  success = hipblasDdot(handle, n, x1_, incx_, y, incy, &resI);
  hipDeviceSynchronize();
  if (HIPBLAS_STATUS_SUCCESS != success)
    throw std::runtime_error("Error: hipblas_dot returned error code.");
  return std::complex<double>{resR, resI};
}

inline hipblasStatus_t hipblas_axpy(hipblasHandle_t handle,
                                    int n,
                                    const float alpha,
                                    const float* x,
                                    int incx,
                                    float* y,
                                    int incy)
{
  hipblasStatus_t success = hipblasSaxpy(handle, n, &alpha, x, incx, y, incy);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_axpy(hipblasHandle_t handle,
                                    int n,
                                    const double alpha,
                                    const double* x,
                                    int incx,
                                    double* y,
                                    int incy)
{
  hipblasStatus_t success = hipblasDaxpy(handle, n, &alpha, x, incx, y, incy);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_axpy(hipblasHandle_t handle,
                                    int n,
                                    const std::complex<float> alpha,
                                    const std::complex<float>* x,
                                    int incx,
                                    std::complex<float>* y,
                                    int incy)
{
  hipblasStatus_t success =
      hipblasCaxpy(handle, n, reinterpret_cast<hipblasComplex const*>(&alpha),
                   reinterpret_cast<hipblasComplex const*>(x), incx, reinterpret_cast<hipblasComplex*>(y), incy);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_axpy(hipblasHandle_t handle,
                                    int n,
                                    const std::complex<double> alpha,
                                    const std::complex<double>* x,
                                    int incx,
                                    std::complex<double>* y,
                                    int incy)
{
  hipblasStatus_t success = hipblasZaxpy(handle, n, reinterpret_cast<hipblasDoubleComplex const*>(&alpha),
                                         reinterpret_cast<hipblasDoubleComplex const*>(x), incx,
                                         reinterpret_cast<hipblasDoubleComplex*>(y), incy);
  hipDeviceSynchronize();
  return success;
}

// Level-2
inline hipblasStatus_t hipblas_gemv(hipblasHandle_t handle,
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
  hipblasStatus_t success =
      hipblasSgemv(handle, hipblasOperation(Atrans), M, N, &alpha, A, lda, x, incx, &beta, y, incy);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_gemv(hipblasHandle_t handle,
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
  hipblasStatus_t success =
      hipblasDgemv(handle, hipblasOperation(Atrans), M, N, &alpha, A, lda, x, incx, &beta, y, incy);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_gemv(hipblasHandle_t handle,
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
  hipblasStatus_t success =
      hipblasCgemv(handle, hipblasOperation(Atrans), M, N, reinterpret_cast<hipblasComplex const*>(&alpha),
                   reinterpret_cast<hipblasComplex const*>(A), lda, reinterpret_cast<hipblasComplex const*>(x), incx,
                   reinterpret_cast<hipblasComplex const*>(&beta), reinterpret_cast<hipblasComplex*>(y), incy);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_gemv(hipblasHandle_t handle,
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
  hipblasStatus_t success =
      hipblasZgemv(handle, hipblasOperation(Atrans), M, N, reinterpret_cast<hipblasDoubleComplex const*>(&alpha),
                   reinterpret_cast<hipblasDoubleComplex const*>(A), lda,
                   reinterpret_cast<hipblasDoubleComplex const*>(x), incx,
                   reinterpret_cast<hipblasDoubleComplex const*>(&beta), reinterpret_cast<hipblasDoubleComplex*>(y),
                   incy);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_gemv(hipblasHandle_t handle,
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
  hipblasStatus_t success;
  char Nt('N');
  char Tt('T');
  if (Atrans == 'n' || Atrans == 'N')
    success =
        hipblasSgemm(handle, hipblasOperation(Nt), hipblasOperation(Tt), 2, M, N, &alpha,
                     reinterpret_cast<float const*>(x), 2 * incx, A, lda, &beta, reinterpret_cast<float*>(y), 2 * incy);
  else if (Atrans == 't' || Atrans == 'T')
    success =
        hipblasSgemm(handle, hipblasOperation(Nt), hipblasOperation(Nt), 2, N, M, &alpha,
                     reinterpret_cast<float const*>(x), 2 * incx, A, lda, &beta, reinterpret_cast<float*>(y), 2 * incy);
  else
    assert(0);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_gemv(hipblasHandle_t handle,
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
  hipblasStatus_t success;
  char Nt('N');
  char Tt('T');
  if (Atrans == 'n' || Atrans == 'N')
    success = hipblasDgemm(handle, hipblasOperation(Nt), hipblasOperation(Tt), 2, M, N, &alpha,
                           reinterpret_cast<double const*>(x), 2 * incx, A, lda, &beta, reinterpret_cast<double*>(y),
                           2 * incy);
  else if (Atrans == 't' || Atrans == 'T')
    success = hipblasDgemm(handle, hipblasOperation(Nt), hipblasOperation(Nt), 2, N, M, &alpha,
                           reinterpret_cast<double const*>(x), 2 * incx, A, lda, &beta, reinterpret_cast<double*>(y),
                           2 * incy);
  else
    assert(0);
  hipDeviceSynchronize();
  return success;
}


// Level-3
inline hipblasStatus_t hipblas_gemm(hipblasHandle_t handle,
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
  hipblasStatus_t success = hipblasSgemm(handle, hipblasOperation(Atrans), hipblasOperation(Btrans), M, N, K, &alpha, A,
                                         lda, B, ldb, &beta, C, ldc);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_gemm(hipblasHandle_t handle,
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
  hipblasStatus_t success = hipblasDgemm(handle, hipblasOperation(Atrans), hipblasOperation(Btrans), M, N, K, &alpha, A,
                                         lda, B, ldb, &beta, C, ldc);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_gemm(hipblasHandle_t handle,
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
  hipblasStatus_t success =
      hipblasCgemm(handle, hipblasOperation(Atrans), hipblasOperation(Btrans), M, N, K,
                   reinterpret_cast<hipblasComplex const*>(&alpha), reinterpret_cast<hipblasComplex const*>(A), lda,
                   reinterpret_cast<hipblasComplex const*>(B), ldb, reinterpret_cast<hipblasComplex const*>(&beta),
                   reinterpret_cast<hipblasComplex*>(C), ldc);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_gemm(hipblasHandle_t handle,
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
  hipblasStatus_t success = hipblasZgemm(handle, hipblasOperation(Atrans), hipblasOperation(Btrans), M, N, K,
                                         reinterpret_cast<hipblasDoubleComplex const*>(&alpha),
                                         reinterpret_cast<hipblasDoubleComplex const*>(A), lda,
                                         reinterpret_cast<hipblasDoubleComplex const*>(B), ldb,
                                         reinterpret_cast<hipblasDoubleComplex const*>(&beta),
                                         reinterpret_cast<hipblasDoubleComplex*>(C), ldc);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_gemm(hipblasHandle_t handle,
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
  hipblasStatus_t success =
      hipblasSgemm(handle, hipblasOperation(Atrans), hipblasOperation(Btrans), 2 * M, N, K, &alpha,
                   reinterpret_cast<float const*>(A), 2 * lda, B, ldb, &beta, reinterpret_cast<float*>(C), 2 * ldc);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_gemm(hipblasHandle_t handle,
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
  hipblasStatus_t success =
      hipblasDgemm(handle, hipblasOperation(Atrans), hipblasOperation(Btrans), 2 * M, N, K, &alpha,
                   reinterpret_cast<double const*>(A), 2 * lda, B, ldb, &beta, reinterpret_cast<double*>(C), 2 * ldc);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_gemm(hipblasHandle_t handle,
                                    char Atrans,
                                    char Btrans,
                                    int M,
                                    int N,
                                    int K,
                                    const hipblasDoubleComplex alpha,
                                    const hipblasDoubleComplex* A,
                                    int lda,
                                    const hipblasDoubleComplex* B,
                                    int ldb,
                                    const hipblasDoubleComplex beta,
                                    hipblasDoubleComplex* C,
                                    int ldc)
{
  hipblasStatus_t success = hipblasZgemm(handle, hipblasOperation(Atrans), hipblasOperation(Btrans), M, N, K, &alpha, A,
                                         lda, B, ldb, &beta, C, ldc);
  hipDeviceSynchronize();
  return success;
}

// Extensions
inline hipblasStatus_t hipblas_getrfBatched(hipblasHandle_t handle,
                                            int n,
                                            float** Aarray,
                                            int lda,
                                            int* PivotArray,
                                            int* infoArray,
                                            int batchSize)
{
  hipblasStatus_t success = hipblasSgetrfBatched(handle, n, Aarray, lda, PivotArray, infoArray, batchSize);
  //hipblasStatus_t success =
  hipDeviceSynchronize();
  return success;
}

// TODO: Update to hipblas call when this gets merged upstream.
inline hipblasStatus_t hipblas_getrfBatched(hipblasHandle_t handle,
                                            int n,
                                            double** Aarray,
                                            int lda,
                                            int* PivotArray,
                                            int ldc,
                                            int* infoArray,
                                            int batchSize)
{
  hipblasStatus_t success = hipblasDgetrfBatched(handle, n, Aarray, lda, PivotArray, infoArray, batchSize);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_getrfBatched(hipblasHandle_t handle,
                                            int n,
                                            std::complex<double>** Aarray,
                                            int lda,
                                            int* PivotArray,
                                            int* infoArray,
                                            int batchSize)
{
  hipblasStatus_t success = hipblasZgetrfBatched(handle, n, reinterpret_cast<hipblasDoubleComplex* const*>(Aarray), lda,
                                                 PivotArray, infoArray, batchSize);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_getrfBatched(hipblasHandle_t handle,
                                            int n,
                                            std::complex<float>** Aarray,
                                            int lda,
                                            int* PivotArray,
                                            int* infoArray,
                                            int batchSize)
{
  hipblasStatus_t success = hipblasCgetrfBatched(handle, n, reinterpret_cast<hipblasComplex* const*>(Aarray), lda,
                                                 PivotArray, infoArray, batchSize);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_getrsBatched(hipblasHandle_t handle,
                                            hipblasOperation_t op,
                                            int n,
                                            int nrhs,
                                            float** Aarray,
                                            int lda,
                                            int* PivotArray,
                                            float** Carray,
                                            int ldc,
                                            int* infoArray,
                                            int batchSize)
{
  hipblasStatus_t success =
      hipblasSgetrsBatched(handle, op, n, nrhs, Aarray, lda, PivotArray, Carray, ldc, infoArray, batchSize);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_getrsBatched(hipblasHandle_t handle,
                                            hipblasOperation_t op,
                                            int n,
                                            int nrhs,
                                            double** Aarray,
                                            int lda,
                                            int* PivotArray,
                                            double** Carray,
                                            int ldc,
                                            int* infoArray,
                                            int batchSize)
{
  hipblasStatus_t success =
      hipblasDgetrsBatched(handle, op, n, nrhs, Aarray, lda, PivotArray, Carray, ldc, infoArray, batchSize);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_getriBatched(hipblasHandle_t handle,
                                            hipblasOperation_t op,
                                            int n,
                                            int nrhs,
                                            std::complex<double>** Aarray,
                                            int lda,
                                            int* PivotArray,
                                            std::complex<double>** Carray,
                                            int ldc,
                                            int* infoArray,
                                            int batchSize)
{
  hipblasStatus_t success =
      hipblasZgetrsBatched(handle, op, n, n, reinterpret_cast<hipblasDoubleComplex* const*>(Aarray), lda, PivotArray,
                           reinterpret_cast<hipblasDoubleComplex* const*>(Carray), ldc, infoArray, batchSize);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_getriBatched(hipblasHandle_t handle,
                                            hipblasOperation_t op,
                                            int n,
                                            int nrhs,
                                            std::complex<float>** Aarray,
                                            int lda,
                                            int* PivotArray,
                                            std::complex<float>** Carray,
                                            int ldc,
                                            int* infoArray,
                                            int batchSize)
{
  hipblasStatus_t success =
      hipblasCgetrsBatched(handle, op, n, n, reinterpret_cast<hipblasComplex* const*>(Aarray), lda, PivotArray,
                           reinterpret_cast<hipblasComplex* const*>(Carray), ldc, infoArray, batchSize);
  hipDeviceSynchronize();
  return success;
}

// TODO: Not implemented
inline hipblasStatus_t hipblas_matinvBatched(hipblasHandle_t handle,
                                             int n,
                                             float** Aarray,
                                             int lda,
                                             float** Carray,
                                             int ldc,
                                             int* infoArray,
                                             int batchSize)
{
  //hipblasStatus_t success =
  //hipblasSmatinvBatched(handle,n,Aarray,lda,Carray,ldc,infoArray,batchSize);
  throw std::runtime_error("Error: matinvBatched doesn't exist.");
  hipDeviceSynchronize();
  hipblasStatus_t success;
  return success;
}

inline hipblasStatus_t hipblas_matinvBatched(hipblasHandle_t handle,
                                             int n,
                                             double** Aarray,
                                             int lda,
                                             double** Carray,
                                             int ldc,
                                             int* infoArray,
                                             int batchSize)
{
  //hipblasStatus_t success =
  //hipblasDmatinvBatched(handle,n,Aarray,lda,Carray,ldc,infoArray,batchSize);
  throw std::runtime_error("Error: matinvBatched doesn't exist.");
  hipDeviceSynchronize();
  hipblasStatus_t success;
  return success;
}

inline hipblasStatus_t hipblas_matinvBatched(hipblasHandle_t handle,
                                             int n,
                                             std::complex<float>** Aarray,
                                             int lda,
                                             std::complex<float>** Carray,
                                             int ldc,
                                             int* infoArray,
                                             int batchSize)
{
  //hipblasStatus_t success =
  //hipblasCmatinvBatched(handle,n,
  //reinterpret_cast<const hipblasComplex * const*>(Aarray),lda,
  //reinterpret_cast<hipblasComplex **>(Carray),ldc,
  //infoArray,batchSize);
  throw std::runtime_error("Error: matinvBatched doesn't exist.");
  hipDeviceSynchronize();
  hipblasStatus_t success;
  return success;
}

inline hipblasStatus_t hipblas_matinvBatched(hipblasHandle_t handle,
                                             int n,
                                             std::complex<double>** Aarray,
                                             int lda,
                                             std::complex<double>** Carray,
                                             int ldc,
                                             int* infoArray,
                                             int batchSize)
{
  //hipblasStatus_t success =
  //hipblasZmatinvBatched(handle,n,
  //reinterpret_cast<const hipblasDoubleComplex *const*>(Aarray),lda,
  //reinterpret_cast<hipblasDoubleComplex **>(Carray),ldc,
  //infoArray,batchSize);
  throw std::runtime_error("Error: matinvBatched doesn't exist.");
  hipDeviceSynchronize();
  hipblasStatus_t success;
  return success;
}

inline hipblasStatus_t hipblas_geam(hipblasHandle_t handle,
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
  hipblasStatus_t success = hipblasSgeam(handle, hipblasOperation(Atrans), hipblasOperation(Btrans), M, N, &alpha, A,
                                         lda, &beta, B, ldb, C, ldc);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_geam(hipblasHandle_t handle,
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
  hipblasStatus_t success = hipblasDgeam(handle, hipblasOperation(Atrans), hipblasOperation(Btrans), M, N, &alpha, A,
                                         lda, &beta, B, ldb, C, ldc);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_geam(hipblasHandle_t handle,
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
  hipblasStatus_t success =
      hipblasCgeam(handle, hipblasOperation(Atrans), hipblasOperation(Btrans), M, N,
                   reinterpret_cast<hipblasComplex const*>(&alpha), reinterpret_cast<hipblasComplex const*>(A), lda,
                   reinterpret_cast<hipblasComplex const*>(&beta), reinterpret_cast<hipblasComplex const*>(B), ldb,
                   reinterpret_cast<hipblasComplex*>(C), ldc);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_geam(hipblasHandle_t handle,
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
  hipblasStatus_t success = hipblasZgeam(handle, hipblasOperation(Atrans), hipblasOperation(Btrans), M, N,
                                         reinterpret_cast<hipblasDoubleComplex const*>(&alpha),
                                         reinterpret_cast<hipblasDoubleComplex const*>(A), lda,
                                         reinterpret_cast<hipblasDoubleComplex const*>(&beta),
                                         reinterpret_cast<hipblasDoubleComplex const*>(B), ldb,
                                         reinterpret_cast<hipblasDoubleComplex*>(C), ldc);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_gemmStridedBatched(hipblasHandle_t handle,
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
  hipblasStatus_t success =
      hipblasSgemmStridedBatched(handle, hipblasOperation(Atrans), hipblasOperation(Btrans), M, N, K, &alpha, A, lda,
                                 strideA, B, ldb, strideB, &beta, C, ldc, strideC, batchSize);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_gemmStridedBatched(hipblasHandle_t handle,
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
  hipblasStatus_t success =
      hipblasDgemmStridedBatched(handle, hipblasOperation(Atrans), hipblasOperation(Btrans), M, N, K, &alpha, A, lda,
                                 strideA, B, ldb, strideB, &beta, C, ldc, strideC, batchSize);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_gemmStridedBatched(hipblasHandle_t handle,
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
  hipblasStatus_t success = hipblasCgemmStridedBatched(handle, hipblasOperation(Atrans), hipblasOperation(Btrans), M, N,
                                                       K, reinterpret_cast<hipblasComplex const*>(&alpha),
                                                       reinterpret_cast<hipblasComplex const*>(A), lda, strideA,
                                                       reinterpret_cast<hipblasComplex const*>(B), ldb, strideB,
                                                       reinterpret_cast<hipblasComplex const*>(&beta),
                                                       reinterpret_cast<hipblasComplex*>(C), ldc, strideC, batchSize);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_gemmStridedBatched(hipblasHandle_t handle,
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
  hipblasStatus_t success =
      hipblasZgemmStridedBatched(handle, hipblasOperation(Atrans), hipblasOperation(Btrans), M, N, K,
                                 reinterpret_cast<hipblasDoubleComplex const*>(&alpha),
                                 reinterpret_cast<hipblasDoubleComplex const*>(A), lda, strideA,
                                 reinterpret_cast<hipblasDoubleComplex const*>(B), ldb, strideB,
                                 reinterpret_cast<hipblasDoubleComplex const*>(&beta),
                                 reinterpret_cast<hipblasDoubleComplex*>(C), ldc, strideC, batchSize);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_gemmBatched(hipblasHandle_t handle,
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
  hipblasStatus_t success = hipblasSgemmBatched(handle, hipblasOperation(Atrans), hipblasOperation(Btrans), M, N, K,
                                                &alpha, A, lda, B, ldb, &beta, C, ldc, batchSize);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_gemmBatched(hipblasHandle_t handle,
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
  hipblasStatus_t success = hipblasDgemmBatched(handle, hipblasOperation(Atrans), hipblasOperation(Btrans), M, N, K,
                                                &alpha, A, lda, B, ldb, &beta, C, ldc, batchSize);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_gemmBatched(hipblasHandle_t handle,
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
  hipblasStatus_t success =
      hipblasCgemmBatched(handle, hipblasOperation(Atrans), hipblasOperation(Btrans), M, N, K,
                          reinterpret_cast<hipblasComplex*>(&alpha), reinterpret_cast<hipblasComplex**>(A), lda,
                          reinterpret_cast<hipblasComplex**>(B), ldb, reinterpret_cast<hipblasComplex*>(&beta),
                          reinterpret_cast<hipblasComplex**>(C), ldc, batchSize);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_gemmBatched(hipblasHandle_t handle,
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
  hipblasStatus_t success =
      hipblasZgemmBatched(handle, hipblasOperation(Atrans), hipblasOperation(Btrans), M, N, K,
                          reinterpret_cast<hipblasDoubleComplex*>(&alpha), reinterpret_cast<hipblasDoubleComplex**>(A),
                          lda, reinterpret_cast<hipblasDoubleComplex**>(B), ldb,
                          reinterpret_cast<hipblasDoubleComplex*>(&beta), reinterpret_cast<hipblasDoubleComplex**>(C),
                          ldc, batchSize);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_gemmBatched(hipblasHandle_t handle,
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
  hipblasStatus_t success = hipblasSgemmBatched(handle, hipblasOperation(Atrans), hipblasOperation(Btrans), 2 * M, N, K,
                                                &alpha, reinterpret_cast<float**>(A), 2 * lda, B, ldb, &beta,
                                                reinterpret_cast<float**>(C), 2 * ldc, batchSize);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_gemmBatched(hipblasHandle_t handle,
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
  hipblasStatus_t success = hipblasDgemmBatched(handle, hipblasOperation(Atrans), hipblasOperation(Btrans), 2 * M, N, K,
                                                &alpha, reinterpret_cast<double**>(A), 2 * lda, B, ldb, &beta,
                                                reinterpret_cast<double**>(C), 2 * ldc, batchSize);
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_geqrfBatched(hipblasHandle_t handle,
                                            int m,
                                            int n,
                                            double** Aarray,
                                            int lda,
                                            double** TauArray,
                                            int* info,
                                            int batchSize)
{
  //hipblasStatus_t success = hipblasDgeqrfBatched(handle,m,n,Aarray,lda,TauArray,info,batchSize);
  hipblasStatus_t success;
  throw std::runtime_error("Error: geqrfBatched doesn't exist.");
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_geqrfBatched(hipblasHandle_t handle,
                                            int m,
                                            int n,
                                            float** Aarray,
                                            int lda,
                                            float** TauArray,
                                            int* info,
                                            int batchSize)
{
  //hipblasStatus_t success = hipblasSgeqrfBatched(handle,m,n,Aarray,lda,TauArray,info,batchSize);
  hipblasStatus_t success;
  throw std::runtime_error("Error: geqrfBatched doesn't exist.");
  hipDeviceSynchronize();
  return success;
}


inline hipblasStatus_t hipblas_geqrfBatched(hipblasHandle_t handle,
                                            int m,
                                            int n,
                                            std::complex<double>** Aarray,
                                            int lda,
                                            std::complex<double>** TauArray,
                                            int* info,
                                            int batchSize)
{
  //hipblasStatus_t success = hipblasZgeqrfBatched(handle,m,n,
  //reinterpret_cast<hipDoubleComplex **>(Aarray),lda,
  //reinterpret_cast<hipDoubleComplex **>(TauArray),info,batchSize);
  hipblasStatus_t success;
  throw std::runtime_error("Error: geqrfBatched doesn't exist.");
  hipDeviceSynchronize();
  return success;
}

inline hipblasStatus_t hipblas_geqrfBatched(hipblasHandle_t handle,
                                            int m,
                                            int n,
                                            std::complex<float>** Aarray,
                                            int lda,
                                            std::complex<float>** TauArray,
                                            int* info,
                                            int batchSize)
{
  //hipblasStatus_t success = hipblasCgeqrfBatched(handle,m,n,
  //reinterpret_cast<hipComplex **>(Aarray),lda,
  //reinterpret_cast<hipComplex **>(TauArray),info,batchSize);
  hipblasStatus_t success;
  throw std::runtime_error("Error: geqrfBatched doesn't exist.");
  hipDeviceSynchronize();
  return success;
}

} // namespace hipblas

#endif
