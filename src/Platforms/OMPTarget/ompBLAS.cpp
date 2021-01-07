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
#include <stdexcept>
#include "config.h"
#if !defined(OPENMP_NO_COMPLEX)
#include "ompReduction.hpp"
#endif

namespace qmcplusplus
{
namespace ompBLAS
{

template<typename T>
ompBLAS_status gemv_impl(ompBLAS_handle& handle,
                         const char      trans,
                         const int       m,
                         const int       n,
                         const T         alpha,
                         const T* const  A,
                         const int       lda,
                         const T* const  x,
                         const int       incx,
                         const T         beta,
                         T* const        y,
                         const int       incy)
{
  if (trans == 'T')
  {
    if (incx !=1 || incy != 1)
      throw std::runtime_error("incx !=1 or incy != 1 are not implemented in ompBLAS::gemv_impl!");

    //BLAS::gemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy);
    PRAGMA_OFFLOAD("omp target teams distribute num_teams(m) is_device_ptr(A, x, y)")
    for(size_t i = 0; i < m; i++)
    {
      T dot_sum(0);
      PRAGMA_OFFLOAD("omp parallel for simd reduction(+: dot_sum)")
      for(size_t j = 0; j < n; j++)
        dot_sum += x[j] * A[i * lda + j];
      y[i] = alpha * dot_sum + beta * y[i];
    }
    return 0;
  }
  else
  {
    throw std::runtime_error("trans = 'N' not implemented in ompBLAS::gemv_impl!");
  }
}

ompBLAS_status gemv(ompBLAS_handle&    handle,
                    const char         trans,
                    const int          m,
                    const int          n,
                    const float        alpha,
                    const float* const A,
                    const int          lda,
                    const float* const x,
                    const int          incx,
                    const float        beta,
                    float* const       y,
                    const int          incy)
{
  return gemv_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy);
}

ompBLAS_status gemv(ompBLAS_handle&     handle,
                    const char          trans,
                    const int           m,
                    const int           n,
                    const double        alpha,
                    const double* const A,
                    const int           lda,
                    const double* const x,
                    const int           incx,
                    const double        beta,
                    double* const       y,
                    const int           incy)
{
  return gemv_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy);
}

#if !defined(OPENMP_NO_COMPLEX)
ompBLAS_status gemv(ompBLAS_handle&                  handle,
                    const char                       trans,
                    const int                        m,
                    const int                        n,
                    const std::complex<float>        alpha,
                    const std::complex<float>* const A,
                    const int                        lda,
                    const std::complex<float>* const x,
                    const int                        incx,
                    const std::complex<float>        beta,
                    std::complex<float>* const       y,
                    const int                        incy)
{
  return gemv_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy);
}

ompBLAS_status gemv(ompBLAS_handle&                   handle,
                    const char                        trans,
                    const int                         m,
                    const int                         n,
                    const std::complex<double>        alpha,
                    const std::complex<double>* const A,
                    const int                         lda,
                    const std::complex<double>* const x,
                    const int                         incx,
                    const std::complex<double>        beta,
                    std::complex<double>* const       y,
                    const int                         incy)
{
  return gemv_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy);
}
#endif


template<typename T>
ompBLAS_status gemv_batched_impl(ompBLAS_handle& handle,
                                 const char      trans,
                                 const int       m,
                                 const int       n,
                                 const T*        alpha,
                                 const T* const  A[],
                                 const int       lda,
                                 const T* const  x[],
                                 const int       incx,
                                 const T*        beta,
                                 T* const        y[],
                                 const int       incy,
                                 const int       batch_count)
{
  if (batch_count == 0) return 0;

  if (trans == 'T')
  {
    if (incx !=1 || incy != 1)
      throw std::runtime_error("incx !=1 or incy != 1 are not implemented in ompBLAS::gemv_batched_impl!");

    PRAGMA_OFFLOAD("omp target teams distribute collapse(2) num_teams(batch_count * m) is_device_ptr(A, x, y, alpha, beta)")
    for(size_t ib = 0; ib < batch_count; ib++)
      for(size_t i = 0; i < m; i++)
      {
        T dot_sum(0);
        PRAGMA_OFFLOAD("omp parallel for simd reduction(+: dot_sum)")
        for(size_t j = 0; j < n; j++)
          dot_sum += x[ib][j] * A[ib][i * lda + j];
        y[ib][i] = alpha[ib] * dot_sum + beta[ib] * y[ib][i];
      }
    return 0;
  }
  else
  {
    throw std::runtime_error("trans = 'N' not implemented in gemv_batched_impl!");
  }
}

ompBLAS_status gemv_batched(ompBLAS_handle&    handle,
                            const char         trans,
                            const int          m,
                            const int          n,
                            const float*       alpha,
                            const float* const A[],
                            const int          lda,
                            const float* const x[],
                            const int          incx,
                            const float*       beta,
                            float* const       y[],
                            const int          incy,
                            const int          batch_count)
{
  return gemv_batched_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy, batch_count);
}

ompBLAS_status gemv_batched(ompBLAS_handle&     handle,
                            const char          trans,
                            const int           m,
                            const int           n,
                            const double*       alpha,
                            const double* const A[],
                            const int           lda,
                            const double* const x[],
                            const int           incx,
                            const double*       beta,
                            double* const       y[],
                            const int           incy,
                            const int           batch_count)
{
  return gemv_batched_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy, batch_count);
}

#if !defined(OPENMP_NO_COMPLEX)
ompBLAS_status gemv_batched(ompBLAS_handle&                  handle,
                            const char                       trans,
                            const int                        m,
                            const int                        n,
                            const std::complex<float>*       alpha,
                            const std::complex<float>* const A[],
                            const int                        lda,
                            const std::complex<float>* const x[],
                            const int                        incx,
                            const std::complex<float>*       beta,
                            std::complex<float>* const       y[],
                            const int                        incy,
                            const int                        batch_count)
{
  return gemv_batched_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy, batch_count);
}

ompBLAS_status gemv_batched(ompBLAS_handle&                   handle,
                            const char                        trans,
                            const int                         m,
                            const int                         n,
                            const std::complex<double>*       alpha,
                            const std::complex<double>* const A[],
                            const int                         lda,
                            const std::complex<double>* const x[],
                            const int                         incx,
                            const std::complex<double>*       beta,
                            std::complex<double>* const       y[],
                            const int                         incy,
                            const int                         batch_count)
{
  return gemv_batched_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy, batch_count);
}
#endif


template<typename T>
ompBLAS_status ger_impl(ompBLAS_handle& handle,
                        const int       m,
                        const int       n,
                        const T         alpha,
                        const T* const  x,
                        const int       incx,
                        const T* const  y,
                        const int       incy,
                        T* const        A,
                        const int       lda)
{
  if (incx !=1 || incy != 1)
    throw std::runtime_error("incx !=1 or incy != 1 are not implemented in ompBLAS::ger_impl!");

  //BLAS::ger(m, n, alpha, x, incx, y, incy, A, lda);
  PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(2) is_device_ptr(A, x, y)")
  for(size_t i = 0; i < n; i++)
    for(size_t j = 0; j < m; j++)
      A[i * lda + j] += alpha * x[j] * y[i];
  return 0;
}

ompBLAS_status ger(ompBLAS_handle&    handle,
                   const int          m,
                   const int          n,
                   const float        alpha,
                   const float* const x,
                   const int          incx,
                   const float* const y,
                   const int          incy,
                   float* const       A,
                   const int          lda)
{
  return ger_impl(handle, m, n, alpha, x, incx, y, incy, A, lda);
}

ompBLAS_status ger(ompBLAS_handle&     handle,
                   const int           m,
                   const int           n,
                   const double        alpha,
                   const double* const x,
                   const int           incx,
                   const double* const y,
                   const int           incy,
                   double* const       A,
                   const int           lda)
{
  return ger_impl(handle, m, n, alpha, x, incx, y, incy, A, lda);
}

#if !defined(OPENMP_NO_COMPLEX)
ompBLAS_status ger(ompBLAS_handle&                  handle,
                   const int                        m,
                   const int                        n,
                   const std::complex<float>        alpha,
                   const std::complex<float>* const x,
                   const int                        incx,
                   const std::complex<float>* const y,
                   const int                        incy,
                   std::complex<float>* const       A,
                   const int                        lda)
{
  return ger_impl(handle, m, n, alpha, x, incx, y, incy, A, lda);
}

ompBLAS_status ger(ompBLAS_handle&                   handle,
                   const int                         m,
                   const int                         n,
                   const std::complex<double>        alpha,
                   const std::complex<double>* const x,
                   const int                         incx,
                   const std::complex<double>* const y,
                   const int                         incy,
                   std::complex<double>* const       A,
                   const int                         lda)
{
  return ger_impl(handle, m, n, alpha, x, incx, y, incy, A, lda);
}
#endif


template<typename T>
ompBLAS_status ger_batched_impl(ompBLAS_handle& handle,
                                const int       m,
                                const int       n,
                                const T*        alpha,
                                const T* const  x[],
                                const int       incx,
                                const T* const  y[],
                                const int       incy,
                                T* const        A[],
                                const int       lda,
                                const int       batch_count)
{
  if (batch_count == 0) return 0;

  if (incx !=1 || incy != 1)
    throw std::runtime_error("incx !=1 or incy != 1 are not implemented in ompBLAS::ger_batched_impl!");

  PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(3) is_device_ptr(A, x, y, alpha)")
  for(size_t ib = 0; ib < batch_count; ib++)
    for(size_t i = 0; i < n; i++)
      for(size_t j = 0; j < m; j++)
        A[ib][i * lda + j] += alpha[ib] * x[ib][j] * y[ib][i];
  return 0;
}

ompBLAS_status ger_batched(ompBLAS_handle&    handle,
                           const int          m,
                           const int          n,
                           const float*       alpha,
                           const float* const x[],
                           const int          incx,
                           const float* const y[],
                           const int          incy,
                           float* const       A[],
                           const int          lda,
                           const int          batch_count)
{
  return ger_batched_impl(handle, m, n, alpha, x, incx, y, incy, A, lda, batch_count);
}

ompBLAS_status ger_batched(ompBLAS_handle&     handle,
                           const int           m,
                           const int           n,
                           const double*       alpha,
                           const double* const x[],
                           const int           incx,
                           const double* const y[],
                           const int           incy,
                           double* const       A[],
                           const int           lda,
                           const int           batch_count)
{
  return ger_batched_impl(handle, m, n, alpha, x, incx, y, incy, A, lda, batch_count);
}

#if !defined(OPENMP_NO_COMPLEX)
ompBLAS_status ger_batched(ompBLAS_handle&                  handle,
                           const int                        m,
                           const int                        n,
                           const std::complex<float>*       alpha,
                           const std::complex<float>* const x[],
                           const int                        incx,
                           const std::complex<float>* const y[],
                           const int                        incy,
                           std::complex<float>* const       A[],
                           const int                        lda,
                           const int                        batch_count)
{
  return ger_batched_impl(handle, m, n, alpha, x, incx, y, incy, A, lda, batch_count);
}

ompBLAS_status ger_batched(ompBLAS_handle&                   handle,
                           const int                         m,
                           const int                         n,
                           const std::complex<double>*       alpha,
                           const std::complex<double>* const x[],
                           const int                         incx,
                           const std::complex<double>* const y[],
                           const int                         incy,
                           std::complex<double>* const       A[],
                           const int                         lda,
                           const int                         batch_count)
{
  return ger_batched_impl(handle, m, n, alpha, x, incx, y, incy, A, lda, batch_count);
}
#endif
} // namespace ompBLAS
} // namespace qmcplusplus
