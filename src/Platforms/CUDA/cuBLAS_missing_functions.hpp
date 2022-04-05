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


#ifndef QMCPLUSPLUS_CUBLAS_MISSING_FUNCTIONS_H
#define QMCPLUSPLUS_CUBLAS_MISSING_FUNCTIONS_H

#include <complex>
#include "config.h"
#include "CUDAruntime.hpp"

namespace qmcplusplus
{
/** Implement selected batched BLAS1/2 calls using CUDA for different data types S/C/D/Z.
 * cuBLAS_MFs stands for missing functions in cuBLAS.
 * 1) column major just like the BLAS fortran API
 * 2) all the functions are asynchronous
 * 3) all the pointer arguments are expected as device pointers.
 * 4) in batched APIs, alpha and beta are **not** scalars but pointers to array of batch size.
 */
namespace cuBLAS_MFs
{
using cuBLAS_MFs_status = cudaError_t;
using cuBLAS_MFs_handle = cudaStream_t;

// BLAS2
/** Xgemv batched API
 * @param handle handle for asynchronous computation
 * @param trans whether A matrices are transposed
 * @param m number of rows in A
 * @param n number of columns in A
 * @param alpha the factor vector of A
 * @param A device array of device pointers of matrices
 * @param lda leading dimension of A
 * @param x device array of device pointers of vector
 * @param incx increment for the elements of x. It cannot be zero.
 * @param beta the factor vector of vector y
 * @param y device array of device pointers of vector
 * @param incy increment for the elements of y. It cannot be zero.
 * @param batch_count batch size
 */
cuBLAS_MFs_status gemv_batched(cuBLAS_MFs_handle& handle,
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

cuBLAS_MFs_status gemv_batched(cuBLAS_MFs_handle& handle,
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

cuBLAS_MFs_status gemv_batched(cuBLAS_MFs_handle& handle,
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

cuBLAS_MFs_status gemv_batched(cuBLAS_MFs_handle& handle,
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

/** Xger batched API
 * @param handle handle for asynchronous computation
 * @param m number of rows in A
 * @param n number of columns in A
 * @param alpha the factor vector of A
 * @param x device array of device pointers of vector
 * @param incx increment for the elements of x. It cannot be zero.
 * @param y device array of device pointers of vector
 * @param incy increment for the elements of y. It cannot be zero.
 * @param A device array of device pointers of matrices
 * @param lda leading dimension of A
 * @param batch_count batch size
 */
cuBLAS_MFs_status ger_batched(cuBLAS_MFs_handle& handle,
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

cuBLAS_MFs_status ger_batched(cuBLAS_MFs_handle& handle,
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

cuBLAS_MFs_status ger_batched(cuBLAS_MFs_handle& handle,
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

cuBLAS_MFs_status ger_batched(cuBLAS_MFs_handle& handle,
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

// BLAS1
/** Xcopy batched API
 * @param handle handle for asynchronous computation
 * @param n number of elements to be copied
 * @param in device array of device pointers of vector
 * @param incx increment for the elements of in. It cannot be zero.
 * @param out device array of device pointers of vector
 * @param incy increment for the elements of out. It cannot be zero.
 * @param batch_count batch size
 */
cuBLAS_MFs_status copy_batched(cudaStream_t& hstream,
                               const int n,
                               const float* const in[],
                               const int incx,
                               float* const out[],
                               const int incy,
                               const int batch_count);

cuBLAS_MFs_status copy_batched(cudaStream_t& hstream,
                               const int n,
                               const double* const in[],
                               const int incx,
                               double* const out[],
                               const int incy,
                               const int batch_count);

cuBLAS_MFs_status copy_batched(cudaStream_t& hstream,
                               const int n,
                               const std::complex<float>* const in[],
                               const int incx,
                               std::complex<float>* const out[],
                               const int incy,
                               const int batch_count);

cuBLAS_MFs_status copy_batched(cudaStream_t& hstream,
                               const int n,
                               const std::complex<double>* const in[],
                               const int incx,
                               std::complex<double>* const out[],
                               const int incy,
                               const int batch_count);

} // namespace cuBLAS_MFs

} // namespace qmcplusplus
#endif // QMCPLUSPLUS_CUBLAS_INHOUSE_H
