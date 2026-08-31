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


#include "ompBLAS.hpp"
#include <cstdint>
#include <stdexcept>
#include "config.h"
#if !defined(OPENMP_NO_COMPLEX)
#include "ompReductionComplex.hpp"
#endif

using uint32_t = std::uint32_t;

namespace qmcplusplus
{
namespace ompBLAS
{

template<typename T>
ompBLAS_status gemm_impl(ompBLAS_handle& handle,
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
                         const int ldc)
{
  if (M == 0 || N == 0 || K == 0)
    return 0;

  if (transa == 'T' && transb == 'N') //A(ji) * B(jk) -> C(ik)
  {
    PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(2) is_device_ptr(A, B, C)")
    for (size_t m = 0; m < M; m++)
      for (size_t n = 0; n < N; n++)
      {
        T sum(0);
        for (size_t k = 0; k < K; k++)
          sum += A[lda * m + k] * B[ldb * n + k];
        C[n * ldc + m] = alpha * sum + (beta == T(0) ? T(0) : C[n * ldc + m] * beta);
      }
  }
  else if (transa == 'T' && transb == 'T')
  {
    PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(2) is_device_ptr(A, B, C)")
    for (size_t m = 0; m < M; m++)
      for (size_t n = 0; n < N; n++)
      {
        T sum(0);
        for (size_t k = 0; k < K; k++)
          sum += A[lda * m + k] * B[ldb * k + n];
        C[n * ldc + m] = alpha * sum + (beta == T(0) ? T(0) : C[n * ldc + m] * beta);
      }
  }
  else if (transa == 'N' && transb == 'T')
  {
    PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(2) is_device_ptr(A, B, C)")
    for (size_t m = 0; m < M; m++)
      for (size_t n = 0; n < N; n++)
      {
        T sum(0);
        for (size_t k = 0; k < K; k++)
          sum += A[lda * k + m] * B[ldb * k + n];
        C[n * ldc + m] = alpha * sum + (beta == T(0) ? T(0) : C[n * ldc + m] * beta);
      }
  }
  else if (transa == 'N' && transb == 'N')
  {
    PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(2) is_device_ptr(A, B, C)")
    for (size_t n = 0; n < N; n++)
      for (size_t m = 0; m < M; m++)
      {
        T sum(0);
        for (size_t k = 0; k < K; k++)
          sum += A[lda * k + m] * B[ldb * n + k];
        C[n * ldc + m] = alpha * sum + (beta == T(0) ? T(0) : C[n * ldc + m] * beta);
      }
  }
  else
    throw std::runtime_error("Error: trans=='C' not yet implemented for ompBLAS::gemm.");

  return 0;
}

template<>
ompBLAS_status gemm<float>(ompBLAS_handle& handle,
                           const char transa,
                           const char transb,
                           const int M,
                           const int N,
                           const int K,
                           const float& alpha,
                           const float* const A,
                           const int lda,
                           const float* const B,
                           const int ldb,
                           const float& beta,
                           float* const C,
                           const int ldc)
{
  return gemm_impl(handle, transa, transb, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

template<>
ompBLAS_status gemm<double>(ompBLAS_handle& handle,
                            const char transa,
                            const char transb,
                            const int M,
                            const int N,
                            const int K,
                            const double& alpha,
                            const double* const A,
                            const int lda,
                            const double* const B,
                            const int ldb,
                            const double& beta,
                            double* const C,
                            const int ldc)
{
  return gemm_impl(handle, transa, transb, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

#if !defined(OPENMP_NO_COMPLEX)
template<>
ompBLAS_status gemm<std::complex<float>>(ompBLAS_handle& handle,
                                         const char transa,
                                         const char transb,
                                         const int M,
                                         const int N,
                                         const int K,
                                         const std::complex<float>& alpha,
                                         const std::complex<float>* const A,
                                         const int lda,
                                         const std::complex<float>* const B,
                                         const int ldb,
                                         const std::complex<float>& beta,
                                         std::complex<float>* const C,
                                         const int ldc)
{
  return gemm_impl(handle, transa, transb, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

template<>
ompBLAS_status gemm<std::complex<double>>(ompBLAS_handle& handle,
                                          const char transa,
                                          const char transb,
                                          const int M,
                                          const int N,
                                          const int K,
                                          const std::complex<double>& alpha,
                                          const std::complex<double>* const A,
                                          const int lda,
                                          const std::complex<double>* const B,
                                          const int ldb,
                                          const std::complex<double>& beta,
                                          std::complex<double>* const C,
                                          const int ldc)
{
  return gemm_impl(handle, transa, transb, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}
#endif

template<typename T>
ompBLAS_status gemm_batched_impl(ompBLAS_handle& handle,
                                 const char transa,
                                 const char transb,
                                 const int M,
                                 const int N,
                                 const int K,
                                 const T alpha,
                                 const T* const Aarray[],
                                 const int lda,
                                 const T* const Barray[],
                                 const int ldb,
                                 const T beta,
                                 T* const Carray[],
                                 const int ldc,
                                 const int batch_count)
{
  if (M == 0 || N == 0 || K == 0 || batch_count == 0)
    return 0;

  if (transa == 'T' && transb == 'N') //A(ji) * B(jk) -> C(ik)
  {
    PRAGMA_OFFLOAD("omp target teams distribute is_device_ptr(Aarray, Barray, Carray)")
    for (size_t iw = 0; iw < batch_count; iw++)
    {
      auto A = Aarray[iw];
      auto B = Barray[iw];
      auto C = Carray[iw];
      PRAGMA_OFFLOAD("omp parallel for collapse(2)")
      for (size_t m = 0; m < M; m++)
        for (size_t n = 0; n < N; n++)
        {
          T sum(0);
          for (size_t k = 0; k < K; k++)
            sum += A[lda * m + k] * B[ldb * n + k];
          C[n * ldc + m] = alpha * sum + (beta == T(0) ? T(0) : C[n * ldc + m] * beta);
        }
    }
  }
  else if (transa == 'T' && transb == 'T')
  {
    PRAGMA_OFFLOAD("omp target teams distribute is_device_ptr(Aarray, Barray, Carray)")
    for (size_t iw = 0; iw < batch_count; iw++)
    {
      auto A = Aarray[iw];
      auto B = Barray[iw];
      auto C = Carray[iw];
      PRAGMA_OFFLOAD("omp parallel for collapse(2)")
      for (size_t m = 0; m < M; m++)
        for (size_t n = 0; n < N; n++)
        {
          T sum(0);
          for (size_t k = 0; k < K; k++)
            sum += A[lda * m + k] * B[ldb * k + n];
          C[n * ldc + m] = alpha * sum + (beta == T(0) ? T(0) : C[n * ldc + m] * beta);
        }
    }
  }
  else if (transa == 'N' && transb == 'T')
  {
    PRAGMA_OFFLOAD("omp target teams distribute is_device_ptr(Aarray, Barray, Carray)")
    for (size_t iw = 0; iw < batch_count; iw++)
    {
      auto A = Aarray[iw];
      auto B = Barray[iw];
      auto C = Carray[iw];
      PRAGMA_OFFLOAD("omp parallel for collapse(2)")
      for (size_t m = 0; m < M; m++)
        for (size_t n = 0; n < N; n++)
        {
          T sum(0);
          for (size_t k = 0; k < K; k++)
            sum += A[lda * k + m] * B[ldb * k + n];
          C[n * ldc + m] = alpha * sum + (beta == T(0) ? T(0) : C[n * ldc + m] * beta);
        }
    }
  }
  else if (transa == 'N' && transb == 'N')
  {
    PRAGMA_OFFLOAD("omp target teams distribute is_device_ptr(Aarray, Barray, Carray)")
    for (size_t iw = 0; iw < batch_count; iw++)
    {
      auto A = Aarray[iw];
      auto B = Barray[iw];
      auto C = Carray[iw];
      PRAGMA_OFFLOAD("omp parallel for collapse(2)")
      for (size_t n = 0; n < N; n++)
        for (size_t m = 0; m < M; m++)
        {
          T sum(0);
          for (size_t k = 0; k < K; k++)
            sum += A[lda * k + m] * B[ldb * n + k];
          C[n * ldc + m] = alpha * sum + (beta == T(0) ? T(0) : C[n * ldc + m] * beta);
        }
    }
  }
  else
    throw std::runtime_error("Error: trans=='C' not yet implemented for ompBLAS::gemm.");

  return 0;
}

template<>
ompBLAS_status gemm_batched<float>(ompBLAS_handle& handle,
                                   const char transa,
                                   const char transb,
                                   const int M,
                                   const int N,
                                   const int K,
                                   const float& alpha,
                                   const float* const A[],
                                   const int lda,
                                   const float* const B[],
                                   const int ldb,
                                   const float& beta,
                                   float* const C[],
                                   const int ldc,
                                   const int batch_count)
{
  return gemm_batched_impl(handle, transa, transb, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc, batch_count);
}

template<>
ompBLAS_status gemm_batched<double>(ompBLAS_handle& handle,
                                    const char transa,
                                    const char transb,
                                    const int M,
                                    const int N,
                                    const int K,
                                    const double& alpha,
                                    const double* const A[],
                                    const int lda,
                                    const double* const B[],
                                    const int ldb,
                                    const double& beta,
                                    double* const C[],
                                    const int ldc,
                                    const int batch_count)
{
  return gemm_batched_impl(handle, transa, transb, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc, batch_count);
}

#if !defined(OPENMP_NO_COMPLEX)
template<>
ompBLAS_status gemm_batched<std::complex<float>>(ompBLAS_handle& handle,
                                                 const char transa,
                                                 const char transb,
                                                 const int M,
                                                 const int N,
                                                 const int K,
                                                 const std::complex<float>& alpha,
                                                 const std::complex<float>* const A[],
                                                 const int lda,
                                                 const std::complex<float>* const B[],
                                                 const int ldb,
                                                 const std::complex<float>& beta,
                                                 std::complex<float>* const C[],
                                                 const int ldc,
                                                 const int batch_count)
{
  return gemm_batched_impl(handle, transa, transb, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc, batch_count);
}

template<>
ompBLAS_status gemm_batched<std::complex<double>>(ompBLAS_handle& handle,
                                                  const char transa,
                                                  const char transb,
                                                  const int M,
                                                  const int N,
                                                  const int K,
                                                  const std::complex<double>& alpha,
                                                  const std::complex<double>* const A[],
                                                  const int lda,
                                                  const std::complex<double>* const B[],
                                                  const int ldb,
                                                  const std::complex<double>& beta,
                                                  std::complex<double>* const C[],
                                                  const int ldc,
                                                  const int batch_count)
{
  return gemm_batched_impl(handle, transa, transb, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc, batch_count);
}
#endif

template<typename T>
ompBLAS_status gemv_impl(ompBLAS_handle& handle,
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
                         const int incy)
{
  if (m == 0 || n == 0)
    return 0;

  if (trans == 'T')
  {
    if (incx != 1 || incy != 1)
      throw std::runtime_error("incx!=1 or incy!=1 are not implemented in ompBLAS::gemv_impl trans='T'!");

    PRAGMA_OFFLOAD("omp target teams distribute num_teams(n) is_device_ptr(A, x, y)")
    for (uint32_t i = 0; i < n; i++)
    {
      T dot_sum(0);
      PRAGMA_OFFLOAD("omp parallel for simd reduction(+: dot_sum)")
      for (uint32_t j = 0; j < m; j++)
        dot_sum += x[j] * A[i * lda + j];
      if (beta == T(0))
        y[i] = alpha * dot_sum; // protecting NaN from y
      else
        y[i] = alpha * dot_sum + beta * y[i];
    }
    return 0;
  }
  else if (trans == 'N')
  {
    if (incx != 1 || incy != 1)
      throw std::runtime_error("incx !=1 or incy != 1 are not implemented in ompBLAS::gemv_impl trans='N'!");

    PRAGMA_OFFLOAD("omp target teams distribute num_teams(m) is_device_ptr(A, x, y)")
    for (uint32_t i = 0; i < m; i++)
    {
      T dot_sum(0);
      PRAGMA_OFFLOAD("omp parallel for simd reduction(+: dot_sum)")
      for (uint32_t j = 0; j < n; j++)
        dot_sum += x[j] * A[j * lda + i];
      if (beta == T(0))
        y[i] = alpha * dot_sum; // protecting NaN from y
      else
        y[i] = alpha * dot_sum + beta * y[i];
    }
    return 0;
  }
  else
    throw std::runtime_error("Error: trans=='C' not yet implemented for ompBLAS::gemv_impl.");
}

template<>
ompBLAS_status gemv<float>(ompBLAS_handle& handle,
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
                           const int incy)
{
  return gemv_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy);
}

template<>
ompBLAS_status gemv<double>(ompBLAS_handle& handle,
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
                            const int incy)
{
  return gemv_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy);
}

#if !defined(OPENMP_NO_COMPLEX)
template<>
ompBLAS_status gemv<std::complex<float>>(ompBLAS_handle& handle,
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
                                         const int incy)
{
  return gemv_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy);
}

template<>
ompBLAS_status gemv<std::complex<double>>(ompBLAS_handle& handle,
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
                                          const int incy)
{
  return gemv_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy);
}
#endif


template<typename T>
ompBLAS_status gemv_batched_impl(ompBLAS_handle& handle,
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
  if (m == 0 || n == 0 || batch_count == 0)
    return 0;

  if (trans == 'T')
  {
    if (incx != 1)
      throw std::runtime_error("incx!=1 are not implemented in ompBLAS::gemv_batched_impl trans='T'!");

    PRAGMA_OFFLOAD("omp target teams distribute collapse(2) num_teams(batch_count * n) \
                               is_device_ptr(A, x, y, alpha, beta)")
    for (uint32_t ib = 0; ib < batch_count; ib++)
      for (uint32_t i = 0; i < n; i++)
      {
        T dot_sum(0);
        PRAGMA_OFFLOAD("omp parallel for simd reduction(+: dot_sum)")
        for (uint32_t j = 0; j < m; j++)
          dot_sum += x[ib][j] * A[ib][i * lda + j];
        if (beta[ib] == T(0))
          y[ib][i * incy] = alpha[ib] * dot_sum; // protecting NaN from y
        else
          y[ib][i * incy] = alpha[ib] * dot_sum + beta[ib] * y[ib][i * incy];
      }
    return 0;
  }
  else if (trans == 'N')
  {
    if (incx != 1)
      throw std::runtime_error("incx!=1 are not implemented in ompBLAS::gemv_batched_impl trans='N'!");

    PRAGMA_OFFLOAD("omp target teams distribute collapse(2) num_teams(batch_count * n) \
                               is_device_ptr(A, x, y, alpha, beta)")
    for (uint32_t ib = 0; ib < batch_count; ib++)
      for (uint32_t i = 0; i < m; i++)
      {
        T dot_sum(0);
        PRAGMA_OFFLOAD("omp parallel for simd reduction(+: dot_sum)")
        for (uint32_t j = 0; j < n; j++)
          dot_sum += x[ib][j] * A[ib][j * lda + i];
        if (beta[ib] == T(0))
          y[ib][i * incy] = alpha[ib] * dot_sum; // protecting NaN from y
        else
          y[ib][i * incy] = alpha[ib] * dot_sum + beta[ib] * y[ib][i * incy];
      }
    return 0;
  }
  else
    throw std::runtime_error("Error: trans=='C' not yet implemented for ompBLAS::gemv_impl.");
}

template<>
ompBLAS_status gemv_batched<float>(ompBLAS_handle& handle,
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

template<>
ompBLAS_status gemv_batched<double>(ompBLAS_handle& handle,
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

#if !defined(OPENMP_NO_COMPLEX)
template<>
ompBLAS_status gemv_batched<std::complex<float>>(ompBLAS_handle& handle,
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
  return gemv_batched_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy, batch_count);
}

template<>
ompBLAS_status gemv_batched<std::complex<double>>(ompBLAS_handle& handle,
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
  return gemv_batched_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy, batch_count);
}
#endif


template<typename T>
ompBLAS_status ger_impl(ompBLAS_handle& handle,
                        const int m,
                        const int n,
                        const T alpha,
                        const T* const x,
                        const int incx,
                        const T* const y,
                        const int incy,
                        T* const A,
                        const int lda)
{
  if (m == 0 || n == 0)
    return 0;

  if (incx != 1 || incy != 1)
    throw std::runtime_error("incx !=1 or incy != 1 are not implemented in ompBLAS::ger_impl!");

  //BLAS::ger(m, n, alpha, x, incx, y, incy, A, lda);
  PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(2) is_device_ptr(A, x, y)")
  for (uint32_t i = 0; i < n; i++)
    for (uint32_t j = 0; j < m; j++)
      A[i * lda + j] += alpha * x[j] * y[i];
  return 0;
}

template<>
ompBLAS_status ger<float>(ompBLAS_handle& handle,
                          const int m,
                          const int n,
                          const float alpha,
                          const float* const x,
                          const int incx,
                          const float* const y,
                          const int incy,
                          float* const A,
                          const int lda)
{
  return ger_impl(handle, m, n, alpha, x, incx, y, incy, A, lda);
}

template<>
ompBLAS_status ger<double>(ompBLAS_handle& handle,
                           const int m,
                           const int n,
                           const double alpha,
                           const double* const x,
                           const int incx,
                           const double* const y,
                           const int incy,
                           double* const A,
                           const int lda)
{
  return ger_impl(handle, m, n, alpha, x, incx, y, incy, A, lda);
}

#if !defined(OPENMP_NO_COMPLEX)
template<>
ompBLAS_status ger<std::complex<float>>(ompBLAS_handle& handle,
                                        const int m,
                                        const int n,
                                        const std::complex<float> alpha,
                                        const std::complex<float>* const x,
                                        const int incx,
                                        const std::complex<float>* const y,
                                        const int incy,
                                        std::complex<float>* const A,
                                        const int lda)
{
  return ger_impl(handle, m, n, alpha, x, incx, y, incy, A, lda);
}

template<>
ompBLAS_status ger<std::complex<double>>(ompBLAS_handle& handle,
                                         const int m,
                                         const int n,
                                         const std::complex<double> alpha,
                                         const std::complex<double>* const x,
                                         const int incx,
                                         const std::complex<double>* const y,
                                         const int incy,
                                         std::complex<double>* const A,
                                         const int lda)
{
  return ger_impl(handle, m, n, alpha, x, incx, y, incy, A, lda);
}
#endif


template<typename T>
ompBLAS_status ger_batched_impl(ompBLAS_handle& handle,
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
  if (m == 0 || n == 0 || batch_count == 0)
    return 0;


  if (incx != 1)
    throw std::runtime_error("incx!=1 are not implemented in ompBLAS::ger_batched_impl!");

  PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(3) is_device_ptr(A, x, y, alpha)")
  for (uint32_t ib = 0; ib < batch_count; ib++)
    for (uint32_t i = 0; i < n; i++)
      for (uint32_t j = 0; j < m; j++)
        A[ib][i * lda + j] += alpha[ib] * x[ib][j] * y[ib][i * incy];
  return 0;
}

template<>
ompBLAS_status ger_batched<float>(ompBLAS_handle& handle,
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

template<>
ompBLAS_status ger_batched<double>(ompBLAS_handle& handle,
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

#if !defined(OPENMP_NO_COMPLEX)
template<>
ompBLAS_status ger_batched<std::complex<float>>(ompBLAS_handle& handle,
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
  return ger_batched_impl(handle, m, n, alpha, x, incx, y, incy, A, lda, batch_count);
}

template<>
ompBLAS_status ger_batched<std::complex<double>>(ompBLAS_handle& handle,
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
  return ger_batched_impl(handle, m, n, alpha, x, incx, y, incy, A, lda, batch_count);
}
#endif


template<typename T>
ompBLAS_status copy_batched_impl(ompBLAS_handle& handle,
                                 const int n,
                                 const T* const x[],
                                 const int incx,
                                 T* const y[],
                                 const int incy,
                                 const int batch_count)
{
  if (n == 0 || batch_count == 0)
    return 0;

  PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(2) is_device_ptr(x, y)")
  for (uint32_t ib = 0; ib < batch_count; ib++)
    for (uint32_t i = 0; i < n; i++)
      y[ib][i * incy] = x[ib][i * incx];
  return 0;
}

template<>
ompBLAS_status copy_batched<float>(ompBLAS_handle& handle,
                                   const int n,
                                   const float* const x[],
                                   const int incx,
                                   float* const y[],
                                   const int incy,
                                   const int batch_count)
{
  return copy_batched_impl(handle, n, x, incx, y, incy, batch_count);
}

template<>
ompBLAS_status copy_batched<double>(ompBLAS_handle& handle,
                                    const int n,
                                    const double* const x[],
                                    const int incx,
                                    double* const y[],
                                    const int incy,
                                    const int batch_count)
{
  return copy_batched_impl(handle, n, x, incx, y, incy, batch_count);
}

#if !defined(OPENMP_NO_COMPLEX)
template<>
ompBLAS_status copy_batched<std::complex<float>>(ompBLAS_handle& handle,
                                                 const int n,
                                                 const std::complex<float>* const x[],
                                                 const int incx,
                                                 std::complex<float>* const y[],
                                                 const int incy,
                                                 const int batch_count)
{
  return copy_batched_impl(handle, n, x, incx, y, incy, batch_count);
}

template<>
ompBLAS_status copy_batched<std::complex<double>>(ompBLAS_handle& handle,
                                                  const int n,
                                                  const std::complex<double>* const x[],
                                                  const int incx,
                                                  std::complex<double>* const y[],
                                                  const int incy,
                                                  const int batch_count)
{
  return copy_batched_impl(handle, n, x, incx, y, incy, batch_count);
}
#endif

template<typename T>
ompBLAS_status copy_batched_offset_impl(ompBLAS_handle& handle,
                                        const int n,
                                        const T* const x[],
                                        const int x_offset,
                                        const int incx,
                                        T* const y[],
                                        const int y_offset,
                                        const int incy,
                                        const int batch_count)
{
  if (n == 0 || batch_count == 0)
    return 0;

  PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(2) is_device_ptr(x, y)")
  for (uint32_t ib = 0; ib < batch_count; ib++)
    for (uint32_t i = 0; i < n; i++)
      y[ib][y_offset + i * incy] = x[ib][x_offset + i * incx];
  return 0;
}

template<>
ompBLAS_status copy_batched_offset<float>(ompBLAS_handle& handle,
                                          const int n,
                                          const float* const x[],
                                          const int x_offset,
                                          const int incx,
                                          float* const y[],
                                          const int y_offset,
                                          const int incy,
                                          const int batch_count)
{
  return copy_batched_offset_impl(handle, n, x, x_offset, incx, y, y_offset, incy, batch_count);
}

template<>
ompBLAS_status copy_batched_offset<double>(ompBLAS_handle& handle,
                                           const int n,
                                           const double* const x[],
                                           const int x_offset,
                                           const int incx,
                                           double* const y[],
                                           const int y_offset,
                                           const int incy,
                                           const int batch_count)
{
  return copy_batched_offset_impl(handle, n, x, x_offset, incx, y, y_offset, incy, batch_count);
}

#if !defined(OPENMP_NO_COMPLEX)
template<>
ompBLAS_status copy_batched_offset<std::complex<float>>(ompBLAS_handle& handle,
                                                        const int n,
                                                        const std::complex<float>* const x[],
                                                        const int x_offset,
                                                        const int incx,
                                                        std::complex<float>* const y[],
                                                        const int y_offset,
                                                        const int incy,
                                                        const int batch_count)
{
  return copy_batched_offset_impl(handle, n, x, x_offset, incx, y, y_offset, incy, batch_count);
}

template<>
ompBLAS_status copy_batched_offset<std::complex<double>>(ompBLAS_handle& handle,
                                                         const int n,
                                                         const std::complex<double>* const x[],
                                                         const int x_offset,
                                                         const int incx,
                                                         std::complex<double>* const y[],
                                                         const int y_offset,
                                                         const int incy,
                                                         const int batch_count)
{
  return copy_batched_offset_impl(handle, n, x, x_offset, incx, y, y_offset, incy, batch_count);
}
#endif

template<typename T>
ompBLAS_status copy_impl(ompBLAS_handle& handle,
                         const int n,
                         const T* const x,
                         const int incx,
                         T* const y,
                         const int incy)
{
  if (n == 0)
    return 0;
  PRAGMA_OFFLOAD("omp target teams distribute parallel for is_device_ptr(x, y)")
  for (size_t i = 0; i < n; i++)
    y[i * incy] = x[i * incx];
  return 0;
}

template<>
ompBLAS_status copy<float>(ompBLAS_handle& handle,
                           const int n,
                           const float* const x,
                           const int incx,
                           float* const y,
                           const int incy)
{
  return copy_impl(handle, n, x, incx, y, incy);
}

template<>
ompBLAS_status copy<double>(ompBLAS_handle& handle,
                            const int n,
                            const double* const x,
                            const int incx,
                            double* const y,
                            const int incy)
{
  return copy_impl(handle, n, x, incx, y, incy);
}

template<>
ompBLAS_status copy<std::complex<float>>(ompBLAS_handle& handle,
                                         const int n,
                                         const std::complex<float>* const x,
                                         const int incx,
                                         std::complex<float>* const y,
                                         const int incy)
{
  return copy_impl(handle, n, x, incx, y, incy);
}

template<>
ompBLAS_status copy<std::complex<double>>(ompBLAS_handle& handle,
                                          const int n,
                                          const std::complex<double>* const x,
                                          const int incx,
                                          std::complex<double>* const y,
                                          const int incy)
{
  return copy_impl(handle, n, x, incx, y, incy);
}
} // namespace ompBLAS
} // namespace qmcplusplus
