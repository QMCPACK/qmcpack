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

#ifndef CUSPARSE_FUNCTIONDEFS_DEPRECATED_H
#define CUSPARSE_FUNCTIONDEFS_DEPRECATED_H

#include <cassert>
#include <cuda_runtime.h>
#include "cusparse.h"
#include "AFQMC/Memory/CUDA/cuda_utilities.h"

#if CUSPARSE_VER_MAJOR > 10
#error
#endif

namespace cusparse
{
using qmc_cuda::cusparseOperation;

// Level-2
inline cusparseStatus_t cusparse_csrmv(cusparseHandle_t handle,
                                       char Atrans,
                                       int m,
                                       int n,
                                       int nnz,
                                       const double alpha,
                                       const cusparseMatDescr_t& descrA,
                                       const double* csrValA,
                                       const int* csrRowPtrA,
                                       const int* csrColIndA,
                                       const double* x,
                                       const double beta,
                                       double* y)

{
  cusparseStatus_t success = cusparseDcsrmv(handle, cusparseOperation(Atrans), m, n, nnz, &alpha, descrA, csrValA,
                                           csrRowPtrA, csrColIndA, x, &beta, y);
  cudaDeviceSynchronize();
  return success;
}

inline cusparseStatus_t cusparse_csrmv(cusparseHandle_t handle,
                                       char Atrans,
                                       int m,
                                       int n,
                                       int nnz,
                                       const float alpha,
                                       const cusparseMatDescr_t& descrA,
                                       const float* csrValA,
                                       const int* csrRowPtrA,
                                       const int* csrColIndA,
                                       const float* x,
                                       const float beta,
                                       float* y)

{
  cusparseStatus_t success = cusparseScsrmv(handle, cusparseOperation(Atrans), m, n, nnz, &alpha, descrA, csrValA,
                                           csrRowPtrA, csrColIndA, x, &beta, y);
  cudaDeviceSynchronize();
  return success;
}

inline cusparseStatus_t cusparse_csrmv(cusparseHandle_t handle,
                                       char Atrans,
                                       int m,
                                       int n,
                                       int nnz,
                                       const std::complex<double> alpha,
                                       const cusparseMatDescr_t& descrA,
                                       const std::complex<double>* csrValA,
                                       const int* csrRowPtrA,
                                       const int* csrColIndA,
                                       const std::complex<double>* x,
                                       const std::complex<double> beta,
                                       std::complex<double>* y)

{
  cusparseStatus_t success =
      cusparseZcsrmv(handle, cusparseOperation(Atrans), m, n, nnz, reinterpret_cast<cuDoubleComplex const*>(&alpha),
                     descrA, reinterpret_cast<cuDoubleComplex const*>(csrValA), csrRowPtrA, csrColIndA,
                     reinterpret_cast<cuDoubleComplex const*>(x), reinterpret_cast<cuDoubleComplex const*>(&beta),
                     reinterpret_cast<cuDoubleComplex*>(y));
  cudaDeviceSynchronize();
  return success;
}


inline cusparseStatus_t cusparse_csrmv(cusparseHandle_t handle,
                                       char Atrans,
                                       int m,
                                       int n,
                                       int nnz,
                                       const std::complex<float> alpha,
                                       const cusparseMatDescr_t& descrA,
                                       const std::complex<float>* csrValA,
                                       const int* csrRowPtrA,
                                       const int* csrColIndA,
                                       const std::complex<float>* x,
                                       const std::complex<float> beta,
                                       std::complex<float>* y)
{
  cusparseStatus_t success =
      cusparseCcsrmv(handle, cusparseOperation(Atrans), m, n, nnz, reinterpret_cast<cuComplex const*>(&alpha), descrA,
                     reinterpret_cast<cuComplex const*>(csrValA), csrRowPtrA, csrColIndA,
                     reinterpret_cast<cuComplex const*>(x), reinterpret_cast<cuComplex const*>(&beta),
                     reinterpret_cast<cuComplex*>(y));
  cudaDeviceSynchronize();
  return success;
}


inline cusparseStatus_t cusparse_csrmm(cusparseHandle_t handle,
                                       char Atrans,
                                       int m,
                                       int n,
                                       int k,
                                       int nnz,
                                       const double alpha,
                                       const cusparseMatDescr_t& descrA,
                                       const double* csrValA,
                                       const int* csrRowPtrA,
                                       const int* csrColIndA,
                                       const double* B,
                                       const int ldb,
                                       const double beta,
                                       double* C,
                                       const int ldc)

{
  cusparseStatus_t success = cusparseDcsrmm(handle, cusparseOperation(Atrans), m, n, k, nnz, &alpha, descrA, csrValA,
                                           csrRowPtrA, csrColIndA, B, ldb, &beta, C, ldc);
  cudaDeviceSynchronize();
  return success;
}

inline cusparseStatus_t cusparse_csrmm(cusparseHandle_t handle,
                                       char Atrans,
                                       int m,
                                       int n,
                                       int k,
                                       int nnz,
                                       const float alpha,
                                       const cusparseMatDescr_t& descrA,
                                       const float* csrValA,
                                       const int* csrRowPtrA,
                                       const int* csrColIndA,
                                       const float* B,
                                       const int ldb,
                                       const float beta,
                                       float* C,
                                       const int ldc)

{
  cusparseStatus_t success = cusparseScsrmm(handle, cusparseOperation(Atrans), m, n, k, nnz, &alpha, descrA, csrValA,
                                           csrRowPtrA, csrColIndA, B, ldb, &beta, C, ldc);
  cudaDeviceSynchronize();
  return success;
}

inline cusparseStatus_t cusparse_csrmm(cusparseHandle_t handle,
                                       char Atrans,
                                       int m,
                                       int n,
                                       int k,
                                       int nnz,
                                       const std::complex<double> alpha,
                                       const cusparseMatDescr_t& descrA,
                                       const std::complex<double>* csrValA,
                                       const int* csrRowPtrA,
                                       const int* csrColIndA,
                                       const std::complex<double>* B,
                                       const int ldb,
                                       const std::complex<double> beta,
                                       std::complex<double>* C,
                                       const int ldc)

{
  cusparseStatus_t success =
      cusparseZcsrmm(handle, cusparseOperation(Atrans), m, n, k, nnz, reinterpret_cast<cuDoubleComplex const*>(&alpha),
                     descrA, reinterpret_cast<cuDoubleComplex const*>(csrValA), csrRowPtrA, csrColIndA,
                     reinterpret_cast<cuDoubleComplex const*>(B), ldb, reinterpret_cast<cuDoubleComplex const*>(&beta),
                     reinterpret_cast<cuDoubleComplex*>(C), ldc);
  cudaDeviceSynchronize();
  return success;
}

inline cusparseStatus_t cusparse_csrmm(cusparseHandle_t handle,
                                       char Atrans,
                                       int m,
                                       int n,
                                       int k,
                                       int nnz,
                                       const std::complex<float> alpha,
                                       const cusparseMatDescr_t& descrA,
                                       const std::complex<float>* csrValA,
                                       const int* csrRowPtrA,
                                       const int* csrColIndA,
                                       const std::complex<float>* B,
                                       const int ldb,
                                       const std::complex<float> beta,
                                       std::complex<float>* C,
                                       const int ldc)

{
  cusparseStatus_t success =
      cusparseCcsrmm(handle, cusparseOperation(Atrans), m, n, k, nnz, reinterpret_cast<cuComplex const*>(&alpha),
                     descrA, reinterpret_cast<cuComplex const*>(csrValA), csrRowPtrA, csrColIndA,
                     reinterpret_cast<cuComplex const*>(B), ldb, reinterpret_cast<cuComplex const*>(&beta),
                     reinterpret_cast<cuComplex*>(C), ldc);
  cudaDeviceSynchronize();
  return success;
}

inline cusparseStatus_t cusparse_csrmm2(cusparseHandle_t handle,
                                        char Atrans,
                                        char Btrans,
                                        int m,
                                        int n,
                                        int k,
                                        int nnz,
                                        const double alpha,
                                        const cusparseMatDescr_t& descrA,
                                        const double* csrValA,
                                        const int* csrRowPtrA,
                                        const int* csrColIndA,
                                        const double* B,
                                        const int ldb,
                                        const double beta,
                                        double* C,
                                        const int ldc)

{
  cusparseStatus_t success = cusparseDcsrmm2(handle, cusparseOperation(Atrans), cusparseOperation(Btrans), m, n, k, nnz,
                                            &alpha, descrA, csrValA, csrRowPtrA, csrColIndA, B, ldb, &beta, C, ldc);
  cudaDeviceSynchronize();
  return success;
}

inline cusparseStatus_t cusparse_csrmm2(cusparseHandle_t handle,
                                        char Atrans,
                                        char Btrans,
                                        int m,
                                        int n,
                                        int k,
                                        int nnz,
                                        const float alpha,
                                        const cusparseMatDescr_t& descrA,
                                        const float* csrValA,
                                        const int* csrRowPtrA,
                                        const int* csrColIndA,
                                        const float* B,
                                        const int ldb,
                                        const float beta,
                                        float* C,
                                        const int ldc)

{
  cusparseStatus_t success = cusparseScsrmm2(handle, cusparseOperation(Atrans), cusparseOperation(Btrans), m, n, k, nnz,
                                            &alpha, descrA, csrValA, csrRowPtrA, csrColIndA, B, ldb, &beta, C, ldc);
  cudaDeviceSynchronize();
  return success;
}

inline cusparseStatus_t cusparse_csrmm2(cusparseHandle_t handle,
                                        char Atrans,
                                        char Btrans,
                                        int m,
                                        int n,
                                        int k,
                                        int nnz,
                                        const std::complex<double> alpha,
                                        const cusparseMatDescr_t& descrA,
                                        const std::complex<double>* csrValA,
                                        const int* csrRowPtrA,
                                        const int* csrColIndA,
                                        const std::complex<double>* B,
                                        const int ldb,
                                        const std::complex<double> beta,
                                        std::complex<double>* C,
                                        const int ldc)

{
  cusparseStatus_t success =
      cusparseZcsrmm2(handle, cusparseOperation(Atrans), cusparseOperation(Btrans), m, n, k, nnz,
                      reinterpret_cast<cuDoubleComplex const*>(&alpha), descrA,
                      reinterpret_cast<cuDoubleComplex const*>(csrValA), csrRowPtrA, csrColIndA,
                      reinterpret_cast<cuDoubleComplex const*>(B), ldb, reinterpret_cast<cuDoubleComplex const*>(&beta),
                      reinterpret_cast<cuDoubleComplex*>(C), ldc);
  cudaDeviceSynchronize();
  return success;
}

inline cusparseStatus_t cusparse_csrmm2(cusparseHandle_t handle,
                                        char Atrans,
                                        char Btrans,
                                        int m,
                                        int n,
                                        int k,
                                        int nnz,
                                        const std::complex<float> alpha,
                                        const cusparseMatDescr_t& descrA,
                                        const std::complex<float>* csrValA,
                                        const int* csrRowPtrA,
                                        const int* csrColIndA,
                                        const std::complex<float>* B,
                                        const int ldb,
                                        const std::complex<float> beta,
                                        std::complex<float>* C,
                                        const int ldc)

{
  cusparseStatus_t success =
      cusparseCcsrmm2(handle, cusparseOperation(Atrans), cusparseOperation(Btrans), m, n, k, nnz,
                      reinterpret_cast<cuComplex const*>(&alpha), descrA, reinterpret_cast<cuComplex const*>(csrValA),
                      csrRowPtrA, csrColIndA, reinterpret_cast<cuComplex const*>(B), ldb,
                      reinterpret_cast<cuComplex const*>(&beta), reinterpret_cast<cuComplex*>(C), ldc);
  cudaDeviceSynchronize();
  return success;
}

inline cusparseStatus_t cusparse_gemmi(cusparseHandle_t handle,
                                       int m,
                                       int n,
                                       int k,
                                       int nnz,
                                       const double alpha,
                                       const double* A,
                                       const int lda,
                                       const double* cscValB,
                                       const int* cscColPtrB,
                                       const int* cscRowIndB,
                                       const double beta,
                                       double* C,
                                       const int ldc)

{
  cusparseStatus_t success =
      cusparseDgemmi(handle, m, n, k, nnz, &alpha, A, lda, cscValB, cscColPtrB, cscRowIndB, &beta, C, ldc);
  cudaDeviceSynchronize();
  return success;
}

inline cusparseStatus_t cusparse_gemmi(cusparseHandle_t handle,
                                       int m,
                                       int n,
                                       int k,
                                       int nnz,
                                       const float alpha,
                                       const float* A,
                                       const int lda,
                                       const float* cscValB,
                                       const int* cscColPtrB,
                                       const int* cscRowIndB,
                                       const float beta,
                                       float* C,
                                       const int ldc)

{
  cusparseStatus_t success =
      cusparseSgemmi(handle, m, n, k, nnz, &alpha, A, lda, cscValB, cscColPtrB, cscRowIndB, &beta, C, ldc);
  cudaDeviceSynchronize();
  return success;
}

inline cusparseStatus_t cusparse_gemmi(cusparseHandle_t handle,
                                       int m,
                                       int n,
                                       int k,
                                       int nnz,
                                       const std::complex<double> alpha,
                                       const std::complex<double>* A,
                                       const int lda,
                                       const std::complex<double>* cscValB,
                                       const int* cscColPtrB,
                                       const int* cscRowIndB,
                                       const std::complex<double> beta,
                                       std::complex<double>* C,
                                       const int ldc)

{
  cusparseStatus_t success =
      cusparseZgemmi(handle, m, n, k, nnz, reinterpret_cast<cuDoubleComplex const*>(&alpha),
                     reinterpret_cast<cuDoubleComplex const*>(A), lda,
                     reinterpret_cast<cuDoubleComplex const*>(cscValB), cscColPtrB, cscRowIndB,
                     reinterpret_cast<cuDoubleComplex const*>(&beta), reinterpret_cast<cuDoubleComplex*>(C), ldc);
  cudaDeviceSynchronize();
  return success;
}

inline cusparseStatus_t cusparse_gemmi(cusparseHandle_t handle,
                                       int m,
                                       int n,
                                       int k,
                                       int nnz,
                                       const std::complex<float> alpha,
                                       const std::complex<float>* A,
                                       const int lda,
                                       const std::complex<float>* cscValB,
                                       const int* cscColPtrB,
                                       const int* cscRowIndB,
                                       const std::complex<float> beta,
                                       std::complex<float>* C,
                                       const int ldc)

{
  cusparseStatus_t success =
      cusparseCgemmi(handle, m, n, k, nnz, reinterpret_cast<cuComplex const*>(&alpha),
                     reinterpret_cast<cuComplex const*>(A), lda, reinterpret_cast<cuComplex const*>(cscValB),
                     cscColPtrB, cscRowIndB, reinterpret_cast<cuComplex const*>(&beta), reinterpret_cast<cuComplex*>(C),
                     ldc);
  cudaDeviceSynchronize();
  return success;
}

} // namespace cusparse

#endif
