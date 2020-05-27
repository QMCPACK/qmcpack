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
#include <cuda_runtime_api.h>

namespace qmcplusplus
{
/** interface to cuBLAS_MFs calls for different data types S/C/D/Z
 */
namespace cuBLAS_MFs
{
typedef cudaError_t cuBLAS_MFs_status;
typedef cudaStream_t cuBLAS_MFs_handle;

// BLAS2
// Xgemv_batched
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

// Xger_batched
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
// Xcopy_batched
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
