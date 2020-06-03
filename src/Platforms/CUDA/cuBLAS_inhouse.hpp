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


#ifndef QMCPLUSPLUS_CUBLAS_INHOUSE_H
#define QMCPLUSPLUS_CUBLAS_INHOUSE_H

#include <complex>
#include <cuda_runtime_api.h>

namespace qmcplusplus
{
/** interface to cuBLAS_inhouse calls for different data types S/C/D/Z
 */
namespace cuBLAS_inhouse
{
typedef cudaError_t cuBLAS_inhouse_status;
typedef cudaStream_t cuBLAS_inhouse_handle;

cuBLAS_inhouse_status gemv_batched(cuBLAS_inhouse_handle& handle,
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

cuBLAS_inhouse_status gemv_batched(cuBLAS_inhouse_handle& handle,
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

cuBLAS_inhouse_status gemv_batched(cuBLAS_inhouse_handle& handle,
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

cuBLAS_inhouse_status gemv_batched(cuBLAS_inhouse_handle& handle,
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

cuBLAS_inhouse_status ger_batched(cuBLAS_inhouse_handle& handle,
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

cuBLAS_inhouse_status ger_batched(cuBLAS_inhouse_handle& handle,
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

cuBLAS_inhouse_status ger_batched(cuBLAS_inhouse_handle& handle,
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

cuBLAS_inhouse_status ger_batched(cuBLAS_inhouse_handle& handle,
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

} // namespace cuBLAS_inhouse

} // namespace qmcplusplus
#endif // QMCPLUSPLUS_CUBLAS_INHOUSE_H
