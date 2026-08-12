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


#include "rocBLAS_MFs.hpp"
#include <stdexcept>
#include <string>
#include <rocblas/rocblas.h>

namespace qmcplusplus
{
namespace rocBLAS_MFs
{
namespace
{
/** Map a rocblas_status onto the hipblasStatus_t the caller expects.
 *
 * The enumerations agree only at 0, 1, 3 and 6, so a cast would mislabel failures
 * and, past rocBLAS 11, yield a value that is not a hipblasStatus_t at all. Copied
 * from hipblasConvertStatus in hipBLAS's amd_detail/hipblas.cpp, which is not
 * re-exported.
 */
hipblasStatus_t convertStatus(const rocblas_status status)
{
  switch (status)
  {
  // rocBLAS reports device-memory-size query outcomes as non-zero; none is an error
  case rocblas_status_success:
  case rocblas_status_size_unchanged:
  case rocblas_status_size_increased:
    return HIPBLAS_STATUS_SUCCESS;
  case rocblas_status_invalid_handle:
    return HIPBLAS_STATUS_NOT_INITIALIZED;
  case rocblas_status_not_implemented:
    return HIPBLAS_STATUS_NOT_SUPPORTED;
  case rocblas_status_invalid_pointer:
  case rocblas_status_invalid_size:
  case rocblas_status_invalid_value:
    return HIPBLAS_STATUS_INVALID_VALUE;
  case rocblas_status_memory_error:
    return HIPBLAS_STATUS_ALLOC_FAILED;
  case rocblas_status_internal_error:
    return HIPBLAS_STATUS_INTERNAL_ERROR;
  default:
    return HIPBLAS_STATUS_UNKNOWN;
  }
}

#ifdef QMC_ROCBLAS_ALPHA_BETA_ARRAYS
rocblas_operation convertOperation(const char trans)
{
  if (trans == 'N' || trans == 'n')
    return rocblas_operation_none;
  else if (trans == 'T' || trans == 't')
    return rocblas_operation_transpose;
  else if (trans == 'C' || trans == 'c')
    return rocblas_operation_conjugate_transpose;
  else
    throw std::runtime_error(
        "rocBLAS_MFs::convertOperation trans can only be 'N', 'T', 'C', 'n', 't', 'c'. Input value is " +
        std::string(1, trans));
}

/** Puts the handle in alpha/beta array mode and restores it on scope exit.
 *
 * The strides apply only under rocblas_pointer_mode_device, and both are modal
 * state on a handle shared with every other BLAS call on this stream. The rocBLAS
 * header: "Restore to value 0 if no longer applicable to later function calls."
 * Safe to restore right after an async launch; the values are baked into the
 * kernel arguments.
 */
class AlphaBetaArrayGuard
{
public:
  AlphaBetaArrayGuard(rocblas_handle handle) : handle_(handle), status_(rocblas_status_success)
  {
    if (rocblas_status st = rocblas_get_pointer_mode(handle_, &prior_mode_); st != rocblas_status_success)
    {
      status_ = st;
      return;
    }
    if (rocblas_status st = rocblas_set_pointer_mode(handle_, rocblas_pointer_mode_device);
        st != rocblas_status_success)
    {
      status_ = st;
      return;
    }
    mode_set_ = true;
    if (rocblas_status st = rocblas_set_batch_alpha_stride(handle_, 1); st != rocblas_status_success)
    {
      status_ = st;
      return;
    }
    if (rocblas_status st = rocblas_set_batch_beta_stride(handle_, 1); st != rocblas_status_success)
    {
      status_ = st;
      return;
    }
  }

  ~AlphaBetaArrayGuard()
  {
    rocblas_set_batch_alpha_stride(handle_, 0);
    rocblas_set_batch_beta_stride(handle_, 0);
    if (mode_set_)
      rocblas_set_pointer_mode(handle_, prior_mode_);
  }

  AlphaBetaArrayGuard(const AlphaBetaArrayGuard&)            = delete;
  AlphaBetaArrayGuard& operator=(const AlphaBetaArrayGuard&) = delete;

  /// rocblas_status_success when the handle is fully in alpha/beta array mode
  rocblas_status status() const { return status_; }

private:
  rocblas_handle handle_;
  rocblas_pointer_mode prior_mode_ = rocblas_pointer_mode_host;
  bool mode_set_                   = false;
  rocblas_status status_;
};

rocblas_status gemv_batched_dispatch(rocblas_handle handle,
                                     rocblas_operation trans,
                                     rocblas_int m,
                                     rocblas_int n,
                                     const float* alpha,
                                     const float* const A[],
                                     rocblas_int lda,
                                     const float* const x[],
                                     rocblas_int incx,
                                     const float* beta,
                                     float* const y[],
                                     rocblas_int incy,
                                     rocblas_int batch_count)
{
  return rocblas_sgemv_batched(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy, batch_count);
}

rocblas_status gemv_batched_dispatch(rocblas_handle handle,
                                     rocblas_operation trans,
                                     rocblas_int m,
                                     rocblas_int n,
                                     const double* alpha,
                                     const double* const A[],
                                     rocblas_int lda,
                                     const double* const x[],
                                     rocblas_int incx,
                                     const double* beta,
                                     double* const y[],
                                     rocblas_int incy,
                                     rocblas_int batch_count)
{
  return rocblas_dgemv_batched(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy, batch_count);
}

rocblas_status gemv_batched_dispatch(rocblas_handle handle,
                                     rocblas_operation trans,
                                     rocblas_int m,
                                     rocblas_int n,
                                     const rocblas_float_complex* alpha,
                                     const rocblas_float_complex* const A[],
                                     rocblas_int lda,
                                     const rocblas_float_complex* const x[],
                                     rocblas_int incx,
                                     const rocblas_float_complex* beta,
                                     rocblas_float_complex* const y[],
                                     rocblas_int incy,
                                     rocblas_int batch_count)
{
  return rocblas_cgemv_batched(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy, batch_count);
}

rocblas_status gemv_batched_dispatch(rocblas_handle handle,
                                     rocblas_operation trans,
                                     rocblas_int m,
                                     rocblas_int n,
                                     const rocblas_double_complex* alpha,
                                     const rocblas_double_complex* const A[],
                                     rocblas_int lda,
                                     const rocblas_double_complex* const x[],
                                     rocblas_int incx,
                                     const rocblas_double_complex* beta,
                                     rocblas_double_complex* const y[],
                                     rocblas_int incy,
                                     rocblas_int batch_count)
{
  return rocblas_zgemv_batched(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy, batch_count);
}

/// T is the rocBLAS-native element type; the public overloads reinterpret_cast
/// std::complex to it, which is layout compatible.
template<typename T>
hipblasStatus_t gemv_batched_impl(hipblasHandle_t handle,
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
                                  const int batch_count)
{
  // cuBLAS_MFs treats an empty problem as a no-op, not an error
  if (batch_count == 0 || m == 0 || n == 0)
    return HIPBLAS_STATUS_SUCCESS;

  const rocblas_operation trans_roc = convertOperation(trans);
  auto handle_roc                   = reinterpret_cast<rocblas_handle>(handle);

  AlphaBetaArrayGuard guard(handle_roc);
  if (guard.status() != rocblas_status_success)
    return convertStatus(guard.status());

  return convertStatus(
      gemv_batched_dispatch(handle_roc, trans_roc, m, n, alpha, A, lda, x, incx, beta, y, incy, batch_count));
}

/// geru rather than gerc: cuBLAS_MFs::ger_batched computes A += alpha * x * y^T
/// without conjugating y (cuBLAS_missing_functions.cu:256).
rocblas_status ger_batched_dispatch(rocblas_handle handle,
                                    rocblas_int m,
                                    rocblas_int n,
                                    const float* alpha,
                                    const float* const x[],
                                    rocblas_int incx,
                                    const float* const y[],
                                    rocblas_int incy,
                                    float* const A[],
                                    rocblas_int lda,
                                    rocblas_int batch_count)
{
  return rocblas_sger_batched(handle, m, n, alpha, x, incx, y, incy, A, lda, batch_count);
}

rocblas_status ger_batched_dispatch(rocblas_handle handle,
                                    rocblas_int m,
                                    rocblas_int n,
                                    const double* alpha,
                                    const double* const x[],
                                    rocblas_int incx,
                                    const double* const y[],
                                    rocblas_int incy,
                                    double* const A[],
                                    rocblas_int lda,
                                    rocblas_int batch_count)
{
  return rocblas_dger_batched(handle, m, n, alpha, x, incx, y, incy, A, lda, batch_count);
}

rocblas_status ger_batched_dispatch(rocblas_handle handle,
                                    rocblas_int m,
                                    rocblas_int n,
                                    const rocblas_float_complex* alpha,
                                    const rocblas_float_complex* const x[],
                                    rocblas_int incx,
                                    const rocblas_float_complex* const y[],
                                    rocblas_int incy,
                                    rocblas_float_complex* const A[],
                                    rocblas_int lda,
                                    rocblas_int batch_count)
{
  return rocblas_cgeru_batched(handle, m, n, alpha, x, incx, y, incy, A, lda, batch_count);
}

rocblas_status ger_batched_dispatch(rocblas_handle handle,
                                    rocblas_int m,
                                    rocblas_int n,
                                    const rocblas_double_complex* alpha,
                                    const rocblas_double_complex* const x[],
                                    rocblas_int incx,
                                    const rocblas_double_complex* const y[],
                                    rocblas_int incy,
                                    rocblas_double_complex* const A[],
                                    rocblas_int lda,
                                    rocblas_int batch_count)
{
  return rocblas_zgeru_batched(handle, m, n, alpha, x, incx, y, incy, A, lda, batch_count);
}

/// ger has no beta; setting that stride too is harmless and keeps one guard for
/// both entry points.
template<typename T>
hipblasStatus_t ger_batched_impl(hipblasHandle_t handle,
                                 const int m,
                                 const int n,
                                 const T* alpha,
                                 const T* const x[],
                                 const int incx,
                                 const T* const y[],
                                 const int incy,
                                 T* const A[],
                                 const int lda,
                                 const int batch_count)
{
  // cuBLAS_MFs treats an empty problem as a no-op, not an error
  if (batch_count == 0 || m == 0 || n == 0)
    return HIPBLAS_STATUS_SUCCESS;

  auto handle_roc = reinterpret_cast<rocblas_handle>(handle);

  AlphaBetaArrayGuard guard(handle_roc);
  if (guard.status() != rocblas_status_success)
    return convertStatus(guard.status());

  return convertStatus(ger_batched_dispatch(handle_roc, m, n, alpha, x, incx, y, incy, A, lda, batch_count));
}
#endif // QMC_ROCBLAS_ALPHA_BETA_ARRAYS

rocblas_status copy_batched_dispatch(rocblas_handle handle,
                                     rocblas_int n,
                                     const float* const in[],
                                     rocblas_int incx,
                                     float* const out[],
                                     rocblas_int incy,
                                     rocblas_int batch_count)
{
  return rocblas_scopy_batched(handle, n, in, incx, out, incy, batch_count);
}

rocblas_status copy_batched_dispatch(rocblas_handle handle,
                                     rocblas_int n,
                                     const double* const in[],
                                     rocblas_int incx,
                                     double* const out[],
                                     rocblas_int incy,
                                     rocblas_int batch_count)
{
  return rocblas_dcopy_batched(handle, n, in, incx, out, incy, batch_count);
}

rocblas_status copy_batched_dispatch(rocblas_handle handle,
                                     rocblas_int n,
                                     const rocblas_float_complex* const in[],
                                     rocblas_int incx,
                                     rocblas_float_complex* const out[],
                                     rocblas_int incy,
                                     rocblas_int batch_count)
{
  return rocblas_ccopy_batched(handle, n, in, incx, out, incy, batch_count);
}

rocblas_status copy_batched_dispatch(rocblas_handle handle,
                                     rocblas_int n,
                                     const rocblas_double_complex* const in[],
                                     rocblas_int incx,
                                     rocblas_double_complex* const out[],
                                     rocblas_int incy,
                                     rocblas_int batch_count)
{
  return rocblas_zcopy_batched(handle, n, in, incx, out, incy, batch_count);
}

/// No alpha, so no handle state to set and no guard needed.
template<typename T>
hipblasStatus_t copy_batched_impl(hipblasHandle_t handle,
                                  const int n,
                                  const T* const in[],
                                  const int incx,
                                  T* const out[],
                                  const int incy,
                                  const int batch_count)
{
  if (batch_count == 0 || n == 0)
    return HIPBLAS_STATUS_SUCCESS;

  return convertStatus(
      copy_batched_dispatch(reinterpret_cast<rocblas_handle>(handle), n, in, incx, out, incy, batch_count));
}

} // namespace

//------------------------------------------------------------------------------
#ifdef QMC_ROCBLAS_ALPHA_BETA_ARRAYS
hipblasStatus_t gemv_batched(hipblasHandle_t handle,
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
                             const int batch_count)
{
  return gemv_batched_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy, batch_count);
}

hipblasStatus_t gemv_batched(hipblasHandle_t handle,
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
                             const int batch_count)
{
  return gemv_batched_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy, batch_count);
}

hipblasStatus_t gemv_batched(hipblasHandle_t handle,
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
                             const int batch_count)
{
  return gemv_batched_impl(handle, trans, m, n, reinterpret_cast<const rocblas_float_complex*>(alpha),
                           reinterpret_cast<const rocblas_float_complex* const*>(A), lda,
                           reinterpret_cast<const rocblas_float_complex* const*>(x), incx,
                           reinterpret_cast<const rocblas_float_complex*>(beta),
                           reinterpret_cast<rocblas_float_complex* const*>(y), incy, batch_count);
}

hipblasStatus_t gemv_batched(hipblasHandle_t handle,
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
                             const int batch_count)
{
  return gemv_batched_impl(handle, trans, m, n, reinterpret_cast<const rocblas_double_complex*>(alpha),
                           reinterpret_cast<const rocblas_double_complex* const*>(A), lda,
                           reinterpret_cast<const rocblas_double_complex* const*>(x), incx,
                           reinterpret_cast<const rocblas_double_complex*>(beta),
                           reinterpret_cast<rocblas_double_complex* const*>(y), incy, batch_count);
}

//------------------------------------------------------------------------------
hipblasStatus_t ger_batched(hipblasHandle_t handle,
                            const int m,
                            const int n,
                            const float* alpha,
                            const float* const x[],
                            const int incx,
                            const float* const y[],
                            const int incy,
                            float* const A[],
                            const int lda,
                            const int batch_count)
{
  return ger_batched_impl(handle, m, n, alpha, x, incx, y, incy, A, lda, batch_count);
}

hipblasStatus_t ger_batched(hipblasHandle_t handle,
                            const int m,
                            const int n,
                            const double* alpha,
                            const double* const x[],
                            const int incx,
                            const double* const y[],
                            const int incy,
                            double* const A[],
                            const int lda,
                            const int batch_count)
{
  return ger_batched_impl(handle, m, n, alpha, x, incx, y, incy, A, lda, batch_count);
}

hipblasStatus_t ger_batched(hipblasHandle_t handle,
                            const int m,
                            const int n,
                            const std::complex<float>* alpha,
                            const std::complex<float>* const x[],
                            const int incx,
                            const std::complex<float>* const y[],
                            const int incy,
                            std::complex<float>* const A[],
                            const int lda,
                            const int batch_count)
{
  return ger_batched_impl(handle, m, n, reinterpret_cast<const rocblas_float_complex*>(alpha),
                          reinterpret_cast<const rocblas_float_complex* const*>(x), incx,
                          reinterpret_cast<const rocblas_float_complex* const*>(y), incy,
                          reinterpret_cast<rocblas_float_complex* const*>(A), lda, batch_count);
}

hipblasStatus_t ger_batched(hipblasHandle_t handle,
                            const int m,
                            const int n,
                            const std::complex<double>* alpha,
                            const std::complex<double>* const x[],
                            const int incx,
                            const std::complex<double>* const y[],
                            const int incy,
                            std::complex<double>* const A[],
                            const int lda,
                            const int batch_count)
{
  return ger_batched_impl(handle, m, n, reinterpret_cast<const rocblas_double_complex*>(alpha),
                          reinterpret_cast<const rocblas_double_complex* const*>(x), incx,
                          reinterpret_cast<const rocblas_double_complex* const*>(y), incy,
                          reinterpret_cast<rocblas_double_complex* const*>(A), lda, batch_count);
}
#endif // QMC_ROCBLAS_ALPHA_BETA_ARRAYS

//------------------------------------------------------------------------------
hipblasStatus_t copy_batched(hipblasHandle_t handle,
                             const int n,
                             const float* const in[],
                             const int incx,
                             float* const out[],
                             const int incy,
                             const int batch_count)
{
  return copy_batched_impl(handle, n, in, incx, out, incy, batch_count);
}

hipblasStatus_t copy_batched(hipblasHandle_t handle,
                             const int n,
                             const double* const in[],
                             const int incx,
                             double* const out[],
                             const int incy,
                             const int batch_count)
{
  return copy_batched_impl(handle, n, in, incx, out, incy, batch_count);
}

hipblasStatus_t copy_batched(hipblasHandle_t handle,
                             const int n,
                             const std::complex<float>* const in[],
                             const int incx,
                             std::complex<float>* const out[],
                             const int incy,
                             const int batch_count)
{
  return copy_batched_impl(handle, n, reinterpret_cast<const rocblas_float_complex* const*>(in), incx,
                           reinterpret_cast<rocblas_float_complex* const*>(out), incy, batch_count);
}

hipblasStatus_t copy_batched(hipblasHandle_t handle,
                             const int n,
                             const std::complex<double>* const in[],
                             const int incx,
                             std::complex<double>* const out[],
                             const int incy,
                             const int batch_count)
{
  return copy_batched_impl(handle, n, reinterpret_cast<const rocblas_double_complex* const*>(in), incx,
                           reinterpret_cast<rocblas_double_complex* const*>(out), incy, batch_count);
}

} // namespace rocBLAS_MFs
} // namespace qmcplusplus
