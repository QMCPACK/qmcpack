//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_OMPBLAS_H
#define QMCPLUSPLUS_OMPBLAS_H

#include <complex>

namespace qmcplusplus
{
/** Implement selected batched and non-batched BLAS2 calls using OpenMP offload for different data types S/C/D/Z
 * 1) column major like the BLAS fortran API
 * 2) all the functions are synchronous, expected to be changed to asynchronous in the future.
 * 3) all the pointer arguments are expected as device pointers.
 * 4) in batched APIs, alpha and beta are **not** scalars but pointers to array of batch size.
 */
namespace ompBLAS
{

using ompBLAS_status = int;
using ompBLAS_handle = int;

ompBLAS_status gemv(ompBLAS_handle& handle,
                    const char trans,
                    const int m,
                    const int n,
                    const float alpha,
                    const float* const A,
                    const int lda,
                    const float* const x,
                    const int incx,
                    const float beta,
                    float* const y,
                    const int incy);

ompBLAS_status gemv(ompBLAS_handle& handle,
                    const char trans,
                    const int m,
                    const int n,
                    const double alpha,
                    const double* const A,
                    const int lda,
                    const double* const x,
                    const int incx,
                    const double beta,
                    double* const y,
                    const int incy);

ompBLAS_status gemv(ompBLAS_handle& handle,
                    const char trans,
                    const int m,
                    const int n,
                    const std::complex<float> alpha,
                    const std::complex<float>* const A,
                    const int lda,
                    const std::complex<float>* const x,
                    const int incx,
                    const std::complex<float> beta,
                    std::complex<float>* const y,
                    const int incy);

ompBLAS_status gemv(ompBLAS_handle& handle,
                    const char trans,
                    const int m,
                    const int n,
                    const std::complex<double> alpha,
                    const std::complex<double>* const A,
                    const int lda,
                    const std::complex<double>* const x,
                    const int incx,
                    const std::complex<double> beta,
                    std::complex<double>* const y,
                    const int incy);

ompBLAS_status gemv_batched(ompBLAS_handle& handle,
                            const char trans,
                            const int m,
                            const int n,
                            const float* alpha,
                            const float* const A[],
                            const int lda,
                            const float* const x[],
                            const int incx,
                            const float* beta,
                            float* const y[],
                            const int incy,
                            const int batch_count);

ompBLAS_status gemv_batched(ompBLAS_handle& handle,
                            const char trans,
                            const int m,
                            const int n,
                            const double* alpha,
                            const double* const A[],
                            const int lda,
                            const double* const x[],
                            const int incx,
                            const double* beta,
                            double* const y[],
                            const int incy,
                            const int batch_count);

ompBLAS_status gemv_batched(ompBLAS_handle& handle,
                            const char trans,
                            const int m,
                            const int n,
                            const std::complex<float>* alpha,
                            const std::complex<float>* const A[],
                            const int lda,
                            const std::complex<float>* const x[],
                            const int incx,
                            const std::complex<float>* beta,
                            std::complex<float>* const y[],
                            const int incy,
                            const int batch_count);

ompBLAS_status gemv_batched(ompBLAS_handle& handle,
                            const char trans,
                            const int m,
                            const int n,
                            const std::complex<double>* alpha,
                            const std::complex<double>* const A[],
                            const int lda,
                            const std::complex<double>* const x[],
                            const int incx,
                            const std::complex<double>* beta,
                            std::complex<double>* const y[],
                            const int incy,
                            const int batch_count);

ompBLAS_status ger(ompBLAS_handle& handle,
                   const int m,
                   const int n,
                   const float alpha,
                   const float* const x,
                   const int incx,
                   const float* const y,
                   const int incy,
                   float* const A,
                   const int lda);

ompBLAS_status ger(ompBLAS_handle& handle,
                   const int m,
                   const int n,
                   const double alpha,
                   const double* const x,
                   const int incx,
                   const double* const y,
                   const int incy,
                   double* const A,
                   const int lda);

ompBLAS_status ger(ompBLAS_handle& handle,
                   const int m,
                   const int n,
                   const std::complex<float> alpha,
                   const std::complex<float>* const x,
                   const int incx,
                   const std::complex<float>* const y,
                   const int incy,
                   std::complex<float>* const A,
                   const int lda);

ompBLAS_status ger(ompBLAS_handle& handle,
                   const int m,
                   const int n,
                   const std::complex<double> alpha,
                   const std::complex<double>* const x,
                   const int incx,
                   const std::complex<double>* const y,
                   const int incy,
                   std::complex<double>* const A,
                   const int lda);

ompBLAS_status ger_batched(ompBLAS_handle& handle,
                           const int m,
                           const int n,
                           const float* alpha,
                           const float* const x[],
                           const int incx,
                           const float* const y[],
                           const int incy,
                           float* const A[],
                           const int lda,
                           const int batch_count);

ompBLAS_status ger_batched(ompBLAS_handle& handle,
                           const int m,
                           const int n,
                           const double* alpha,
                           const double* const x[],
                           const int incx,
                           const double* const y[],
                           const int incy,
                           double* const A[],
                           const int lda,
                           const int batch_count);

ompBLAS_status ger_batched(ompBLAS_handle& handle,
                           const int m,
                           const int n,
                           const std::complex<float>* alpha,
                           const std::complex<float>* const x[],
                           const int incx,
                           const std::complex<float>* const y[],
                           const int incy,
                           std::complex<float>* const A[],
                           const int lda,
                           const int batch_count);

ompBLAS_status ger_batched(ompBLAS_handle& handle,
                           const int m,
                           const int n,
                           const std::complex<double>* alpha,
                           const std::complex<double>* const x[],
                           const int incx,
                           const std::complex<double>* const y[],
                           const int incy,
                           std::complex<double>* const A[],
                           const int lda,
                           const int batch_count);

} // namespace ompBLAS

} // namespace qmcplusplus
#endif // QMCPLUSPLUS_OMPBLAS_H
