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


#include "rocBLAS.hpp"
#include <stdexcept>
#include <string>

namespace qmcplusplus
{
namespace rocBLAS
{
#ifdef QMC_ROCBLAS_ALPHA_BETA_ARRAYS
namespace
{
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
        "rocBLAS::convertOperation trans can only be 'N', 'T', 'C', 'n', 't', 'c'. Input value is " +
        std::string(1, trans));
}

/** Puts the handle in alpha/beta array mode and restores it on scope exit.
 *
 * The strides apply only under rocblas_pointer_mode_device, and both are modal
 * state on a handle shared with every other BLAS call on this stream. The rocBLAS
 * header: "Restore to value 0 if no longer applicable to later function calls."
 * Safe to restore right after an async launch; the values are baked into the
 * kernel arguments.
 *
 * ger has no beta; setting that stride too is harmless and keeps one guard for
 * both entry points.
 *
 * This mutates state on the handle for the duration of the call, so the same
 * handle must not be used concurrently from another host thread. QMCPACK gives
 * each Queue its own BLASHandle, so nothing shares one today.
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
    alpha_stride_set_ = true;
    if (rocblas_status st = rocblas_set_batch_beta_stride(handle_, 1); st != rocblas_status_success)
    {
      status_ = st;
      return;
    }
    beta_stride_set_ = true;
  }

  /// undoes what the constructor got as far as, in reverse order
  ~AlphaBetaArrayGuard()
  {
    if (beta_stride_set_)
      rocblas_set_batch_beta_stride(handle_, 0);
    if (alpha_stride_set_)
      rocblas_set_batch_alpha_stride(handle_, 0);
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
  bool alpha_stride_set_           = false;
  bool beta_stride_set_            = false;
  rocblas_status status_;
};
} // namespace

//------------------------------------------------------------------------------
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
                            const int batch_count)
{
  // cuBLAS_MFs treats an empty problem as a no-op, not an error
  if (batch_count == 0 || m == 0 || n == 0)
    return rocblas_status_success;

  AlphaBetaArrayGuard guard(handle);
  if (guard.status() != rocblas_status_success)
    return guard.status();

  return rocblas_sgemv_batched(handle, convertOperation(trans), m, n, alpha, A, lda, x, incx, beta, y, incy,
                               batch_count);
}

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
                            const int batch_count)
{
  if (batch_count == 0 || m == 0 || n == 0)
    return rocblas_status_success;

  AlphaBetaArrayGuard guard(handle);
  if (guard.status() != rocblas_status_success)
    return guard.status();

  return rocblas_dgemv_batched(handle, convertOperation(trans), m, n, alpha, A, lda, x, incx, beta, y, incy,
                               batch_count);
}

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
                            const int batch_count)
{
  if (batch_count == 0 || m == 0 || n == 0)
    return rocblas_status_success;

  AlphaBetaArrayGuard guard(handle);
  if (guard.status() != rocblas_status_success)
    return guard.status();

  // std::complex is layout compatible with the rocBLAS complex types
  return rocblas_cgemv_batched(handle, convertOperation(trans), m, n,
                               reinterpret_cast<const rocblas_float_complex*>(alpha),
                               reinterpret_cast<const rocblas_float_complex* const*>(A), lda,
                               reinterpret_cast<const rocblas_float_complex* const*>(x), incx,
                               reinterpret_cast<const rocblas_float_complex*>(beta),
                               reinterpret_cast<rocblas_float_complex* const*>(y), incy, batch_count);
}

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
                            const int batch_count)
{
  if (batch_count == 0 || m == 0 || n == 0)
    return rocblas_status_success;

  AlphaBetaArrayGuard guard(handle);
  if (guard.status() != rocblas_status_success)
    return guard.status();

  return rocblas_zgemv_batched(handle, convertOperation(trans), m, n,
                               reinterpret_cast<const rocblas_double_complex*>(alpha),
                               reinterpret_cast<const rocblas_double_complex* const*>(A), lda,
                               reinterpret_cast<const rocblas_double_complex* const*>(x), incx,
                               reinterpret_cast<const rocblas_double_complex*>(beta),
                               reinterpret_cast<rocblas_double_complex* const*>(y), incy, batch_count);
}

//------------------------------------------------------------------------------
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
                           const int batch_count)
{
  // cuBLAS_MFs treats an empty problem as a no-op, not an error
  if (batch_count == 0 || m == 0 || n == 0)
    return rocblas_status_success;

  AlphaBetaArrayGuard guard(handle);
  if (guard.status() != rocblas_status_success)
    return guard.status();

  return rocblas_sger_batched(handle, m, n, alpha, x, incx, y, incy, A, lda, batch_count);
}

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
                           const int batch_count)
{
  if (batch_count == 0 || m == 0 || n == 0)
    return rocblas_status_success;

  AlphaBetaArrayGuard guard(handle);
  if (guard.status() != rocblas_status_success)
    return guard.status();

  return rocblas_dger_batched(handle, m, n, alpha, x, incx, y, incy, A, lda, batch_count);
}

/// geru rather than gerc: cuBLAS_MFs::ger_batched computes A += alpha * x * y^T
/// without conjugating y (cuBLAS_missing_functions.cu:256).
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
                           const int batch_count)
{
  if (batch_count == 0 || m == 0 || n == 0)
    return rocblas_status_success;

  AlphaBetaArrayGuard guard(handle);
  if (guard.status() != rocblas_status_success)
    return guard.status();

  return rocblas_cgeru_batched(handle, m, n, reinterpret_cast<const rocblas_float_complex*>(alpha),
                               reinterpret_cast<const rocblas_float_complex* const*>(x), incx,
                               reinterpret_cast<const rocblas_float_complex* const*>(y), incy,
                               reinterpret_cast<rocblas_float_complex* const*>(A), lda, batch_count);
}

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
                           const int batch_count)
{
  if (batch_count == 0 || m == 0 || n == 0)
    return rocblas_status_success;

  AlphaBetaArrayGuard guard(handle);
  if (guard.status() != rocblas_status_success)
    return guard.status();

  return rocblas_zgeru_batched(handle, m, n, reinterpret_cast<const rocblas_double_complex*>(alpha),
                               reinterpret_cast<const rocblas_double_complex* const*>(x), incx,
                               reinterpret_cast<const rocblas_double_complex* const*>(y), incy,
                               reinterpret_cast<rocblas_double_complex* const*>(A), lda, batch_count);
}
#endif // QMC_ROCBLAS_ALPHA_BETA_ARRAYS

//------------------------------------------------------------------------------
rocblas_status copy_batched(rocblas_handle handle,
                            const int n,
                            const float* const in[],
                            const int incx,
                            float* const out[],
                            const int incy,
                            const int batch_count)
{
  // cuBLAS_MFs treats an empty problem as a no-op, not an error
  if (batch_count == 0 || n == 0)
    return rocblas_status_success;

  return rocblas_scopy_batched(handle, n, in, incx, out, incy, batch_count);
}

rocblas_status copy_batched(rocblas_handle handle,
                            const int n,
                            const double* const in[],
                            const int incx,
                            double* const out[],
                            const int incy,
                            const int batch_count)
{
  if (batch_count == 0 || n == 0)
    return rocblas_status_success;

  return rocblas_dcopy_batched(handle, n, in, incx, out, incy, batch_count);
}

rocblas_status copy_batched(rocblas_handle handle,
                            const int n,
                            const std::complex<float>* const in[],
                            const int incx,
                            std::complex<float>* const out[],
                            const int incy,
                            const int batch_count)
{
  if (batch_count == 0 || n == 0)
    return rocblas_status_success;

  return rocblas_ccopy_batched(handle, n, reinterpret_cast<const rocblas_float_complex* const*>(in), incx,
                               reinterpret_cast<rocblas_float_complex* const*>(out), incy, batch_count);
}

rocblas_status copy_batched(rocblas_handle handle,
                            const int n,
                            const std::complex<double>* const in[],
                            const int incx,
                            std::complex<double>* const out[],
                            const int incy,
                            const int batch_count)
{
  if (batch_count == 0 || n == 0)
    return rocblas_status_success;

  return rocblas_zcopy_batched(handle, n, reinterpret_cast<const rocblas_double_complex* const*>(in), incx,
                               reinterpret_cast<rocblas_double_complex* const*>(out), incy, batch_count);
}

} // namespace rocBLAS
} // namespace qmcplusplus
