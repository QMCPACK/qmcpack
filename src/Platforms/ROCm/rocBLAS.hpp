//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
// Copyright(C) 2026 Advanced Micro Devices, Inc. All rights reserved.
//
// File developed by: Jakub Kurzak, jakurzak@amd.com, Advanced Micro Devices, Inc.
//
// File created by: Jakub Kurzak, jakurzak@amd.com, Advanced Micro Devices, Inc.
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_ROCBLAS_H
#define QMCPLUSPLUS_ROCBLAS_H

#include <complex>
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <rocblas/rocblas.h>
#include "config.h"

#define rocblasErrorCheck(ans, cause)                \
  {                                                  \
    rocblasAssert((ans), cause, __FILE__, __LINE__); \
  }
/// prints rocBLAS error messages. Always use rocblasErrorCheck macro.
inline void rocblasAssert(rocblas_status code, const std::string& cause, const char* file, int line)
{
  if (code != rocblas_status_success)
  {
    std::ostringstream err;
    err << "rocblasAssert: " << rocblas_status_to_string(code) << ", file " << file << " , line " << line << std::endl
        << cause << std::endl;
    std::cerr << err.str();
    throw std::runtime_error(cause);
  }
}

namespace qmcplusplus
{
/** Batched rocBLAS calls that hipBLAS does not re-export, used in place of the
 * CUDA/cuBLAS_missing_functions.cu kernels. Same conventions as cuBLAS_MFs:
 * 1) column major just like the BLAS fortran API
 * 2) all the functions are asynchronous
 * 3) all the pointer arguments are expected as device pointers.
 * 4) in batched APIs, alpha and beta are **not** scalars but pointers to array of batch size.
 *
 * Point 4) needs rocblas_set_batch_alpha_stride / _beta_stride, hence the
 * QMC_ROCBLAS_ALPHA_BETA_ARRAYS gate on gemv_batched and ger_batched, an opt-in that is
 * OFF by default. copy_batched has no alpha and is always available.
 *
 * The handle is the rocBLAS handle underlying the hipBLAS handle owned by
 * BLASHandle<PlatformKind::CUDA>; the caller casts it. hipBLAS's AMD backend stores the
 * rocBLAS handle verbatim, so the cast is sound and the calls inherit the stream
 * hipblasSetStream already bound. Same idiom as hipBLAS.cpp.
 */
namespace rocBLAS
{
// BLAS2
#ifdef QMC_ROCBLAS_ALPHA_BETA_ARRAYS
/** Xgemv batched API
 * @param handle rocBLAS handle carrying the stream for asynchronous computation
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
rocblas_status gemv_batched(rocblas_handle handle,
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

rocblas_status gemv_batched(rocblas_handle handle,
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

rocblas_status gemv_batched(rocblas_handle handle,
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

rocblas_status gemv_batched(rocblas_handle handle,
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
 * @param handle rocBLAS handle carrying the stream for asynchronous computation
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
 *
 * Complex maps to rocblas_Xgeru_batched, not gerc: the replaced kernel does not
 * conjugate y.
 */
rocblas_status ger_batched(rocblas_handle handle,
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

rocblas_status ger_batched(rocblas_handle handle,
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

rocblas_status ger_batched(rocblas_handle handle,
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

rocblas_status ger_batched(rocblas_handle handle,
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
#endif // QMC_ROCBLAS_ALPHA_BETA_ARRAYS

// BLAS1
/** Xcopy batched API
 * @param handle rocBLAS handle carrying the stream for asynchronous computation
 * @param n number of elements to be copied
 * @param in device array of device pointers of vector
 * @param incx increment for the elements of in. It cannot be zero.
 * @param out device array of device pointers of vector
 * @param incy increment for the elements of out. It cannot be zero.
 * @param batch_count batch size
 */
rocblas_status copy_batched(rocblas_handle handle,
                            const int n,
                            const float* const in[],
                            const int incx,
                            float* const out[],
                            const int incy,
                            const int batch_count);

rocblas_status copy_batched(rocblas_handle handle,
                            const int n,
                            const double* const in[],
                            const int incx,
                            double* const out[],
                            const int incy,
                            const int batch_count);

rocblas_status copy_batched(rocblas_handle handle,
                            const int n,
                            const std::complex<float>* const in[],
                            const int incx,
                            std::complex<float>* const out[],
                            const int incy,
                            const int batch_count);

rocblas_status copy_batched(rocblas_handle handle,
                            const int n,
                            const std::complex<double>* const in[],
                            const int incx,
                            std::complex<double>* const out[],
                            const int incy,
                            const int batch_count);

} // namespace rocBLAS

} // namespace qmcplusplus
#endif // QMCPLUSPLUS_ROCBLAS_H
