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

template<typename T>
ompBLAS_status gemm(ompBLAS_handle& handle,
                    const char transa,
                    const char transb,
                    const int M,
                    const int N,
                    const int K,
                    const T& alpha,
                    const T* const A,
                    const int lda,
                    const T* const B,
                    const int ldb,
                    const T& beta,
                    T* const C,
                    const int ldc);

template<typename T>
ompBLAS_status gemm_batched(ompBLAS_handle& handle,
                            const char transa,
                            const char transb,
                            const int M,
                            const int N,
                            const int K,
                            const T& alpha,
                            const T* const A[],
                            const int lda,
                            const T* const B[],
                            const int ldb,
                            const T& beta,
                            T* const C[],
                            const int ldc,
                            const int batch_count);

template<typename T>
ompBLAS_status gemv(ompBLAS_handle& handle,
                    const char trans,
                    const int m,
                    const int n,
                    const T alpha,
                    const T* const A,
                    const int lda,
                    const T* const x,
                    const int incx,
                    const T beta,
                    T* const y,
                    const int incy);

template<typename T>
ompBLAS_status gemv_batched(ompBLAS_handle& handle,
                            const char trans,
                            const int m,
                            const int n,
                            const T* alpha,
                            const T* const A[],
                            const int lda,
                            const T* const x[],
                            const int incx,
                            const T* beta,
                            T* const y[],
                            const int incy,
                            const int batch_count);

template<typename T>
ompBLAS_status ger(ompBLAS_handle& handle,
                   const int m,
                   const int n,
                   const T alpha,
                   const T* const x,
                   const int incx,
                   const T* const y,
                   const int incy,
                   T* const A,
                   const int lda);

template<typename T>
ompBLAS_status ger_batched(ompBLAS_handle& handle,
                           const int m,
                           const int n,
                           const T* alpha,
                           const T* const x[],
                           const int incx,
                           const T* const y[],
                           const int incy,
                           T* const A[],
                           const int lda,
                           const int batch_count);

/**
 * @brief copy device data from x to y
 *
 * for b_i in [0,batch_count)
 *   for i   in [0,n)
 *     y[b_i][i*incy] = x[b_i][i*incx]
 *
 * @param n number of elements to copy for each group in the batch
 * @param x,y arrays with length `batch_count`; device pointers to start of data to be copied from(x)/to(y)
 * @param incx,incy storage spacing between elements of x/y to be copied from/to
 * @param batch_count number of batches to process
 */
template<typename T>
ompBLAS_status copy_batched(ompBLAS_handle& handle,
                            const int n,
                            const T* const x[],
                            const int incx,
                            T* const y[],
                            const int incy,
                            const int batch_count);

/**
 * @brief copy device data from x to y with additional offset applied to array of device pointers
 *
 * for b_i in [0,batch_count)
 *   for i   in [0,n)
 *     y[b_i][y_offset + i*incy] = x[b_i][x_offset + i*incx]
 *
 * useful for copying from/to a single row/column of a batch of matrices when a list of device pointers
 * to the start of the matrices is already available
 *
 * @param n number of elements to copy for each group in the batch
 * @param x,y arrays with length `batch_count`; device pointers to start of data to be copied from(x)/to(y)
 * @param x_offset,y_offset distance (in number of elements) from pointer given in x/y to location of first element to be copied
 * @param incx,incy storage spacing between elements of x/y to be copied from/to
 * @param batch_count number of batches to process
 */
template<typename T>
ompBLAS_status copy_batched_offset(ompBLAS_handle& handle,
                                   const int n,
                                   const T* const x[],
                                   const int x_offset,
                                   const int incx,
                                   T* const y[],
                                   const int y_offset,
                                   const int incy,
                                   const int batch_count);

template<typename T>
ompBLAS_status copy(ompBLAS_handle& handle, const int n, const T* const x, const int incx, T* const y, const int incy);

} // namespace ompBLAS

} // namespace qmcplusplus
#endif // QMCPLUSPLUS_OMPBLAS_H
