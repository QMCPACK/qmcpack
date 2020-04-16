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


#include "Platforms/OpenMP/ompBLAS.hpp"
#include <stdexcept>
#include "Numerics/OhmmsBlas.h"

namespace qmcplusplus
{
namespace ompBLAS
{

template<typename T>
ompBLAS_status gemv_impl(ompBLAS_handle& handle,
                         const char&     trans,
                         const int&      m,
                         const int&      n,
                         const T&        alpha,
                         const T* const  A,
                         const int&      lda,
                         const T* const  x,
                         const int&      incx,
                         const T&        beta,
                         T* const        y,
                         const int&      incy)
{
  if (trans == 'T')
  {
    BLAS::gemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy);
    return 0;
  }
  else
  {
    throw std::runtime_error("trans = 'N' not implemented in ompBLAS::gemv_impl!");
  }
}

ompBLAS_status gemv(ompBLAS_handle&    handle,
                    const char&        trans,
                    const int&         m,
                    const int&         n,
                    const float&       alpha,
                    const float* const A,
                    const int&         lda,
                    const float* const x,
                    const int&         incx,
                    const float&       beta,
                    float* const       y,
                    const int&         incy)
{
  return gemv_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy);
}

ompBLAS_status gemv(ompBLAS_handle&     handle,
                    const char&         trans,
                    const int&          m,
                    const int&          n,
                    const double&       alpha,
                    const double* const A,
                    const int&          lda,
                    const double* const x,
                    const int&          incx,
                    const double&       beta,
                    double* const       y,
                    const int&          incy)
{
  return gemv_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy);
}

ompBLAS_status gemv(ompBLAS_handle&                  handle,
                    const char&                      trans,
                    const int&                       m,
                    const int&                       n,
                    const std::complex<float>&       alpha,
                    const std::complex<float>* const A,
                    const int&                       lda,
                    const std::complex<float>* const x,
                    const int&                       incx,
                    const std::complex<float>&       beta,
                    std::complex<float>* const       y,
                    const int&                       incy)
{
  return gemv_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy);
}

ompBLAS_status gemv(ompBLAS_handle&                   handle,
                    const char&                       trans,
                    const int&                        m,
                    const int&                        n,
                    const std::complex<double>&       alpha,
                    const std::complex<double>* const A,
                    const int&                        lda,
                    const std::complex<double>* const x,
                    const int&                        incx,
                    const std::complex<double>&       beta,
                    std::complex<double>* const       y,
                    const int&                        incy)
{
  return gemv_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy);
}


template<typename T>
ompBLAS_status gemv_batched_impl(ompBLAS_handle& handle,
                                 const char&     trans,
                                 const int&      m,
                                 const int&      n,
                                 const T&        alpha,
                                 const T* const  A[],
                                 const int&      lda,
                                 const T* const  x[],
                                 const int&      incx,
                                 const T&        beta,
                                 T* const        y[],
                                 const int&      incy,
                                 const int&      batch_count)
{
  if (trans == 'T')
  {
    return 0;
  }
  else
  {
    throw std::runtime_error("trans = 'N' not implemented in gemv_batched_impl!");
  }
}

ompBLAS_status gemv_batched(ompBLAS_handle&    handle,
                            const char&        trans,
                            const int&         m,
                            const int&         n,
                            const float&       alpha,
                            const float* const A[],
                            const int&         lda,
                            const float* const x[],
                            const int&         incx,
                            const float&       beta,
                            float* const       y[],
                            const int&         incy,
                            const int&         batch_count)
{
  return gemv_batched_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy, batch_count);
}

ompBLAS_status gemv_batched(ompBLAS_handle&     handle,
                            const char&         trans,
                            const int&          m,
                            const int&          n,
                            const double&       alpha,
                            const double* const A[],
                            const int&          lda,
                            const double* const x[],
                            const int&          incx,
                            const double&       beta,
                            double* const       y[],
                            const int&          incy,
                            const int&          batch_count)
{
  return gemv_batched_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy, batch_count);
}

ompBLAS_status gemv_batched(ompBLAS_handle&                  handle,
                            const char&                      trans,
                            const int&                       m,
                            const int&                       n,
                            const std::complex<float>&       alpha,
                            const std::complex<float>* const A[],
                            const int&                       lda,
                            const std::complex<float>* const x[],
                            const int&                       incx,
                            const std::complex<float>&       beta,
                            std::complex<float>* const       y[],
                            const int&                       incy,
                            const int&                       batch_count)
{
  return gemv_batched_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy, batch_count);
}

ompBLAS_status gemv_batched(ompBLAS_handle&                   handle,
                            const char&                       trans,
                            const int&                        m,
                            const int&                        n,
                            const std::complex<double>&       alpha,
                            const std::complex<double>* const A[],
                            const int&                        lda,
                            const std::complex<double>* const x[],
                            const int&                        incx,
                            const std::complex<double>&       beta,
                            std::complex<double>* const       y[],
                            const int&                        incy,
                            const int&                        batch_count)
{
  return gemv_batched_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy, batch_count);
}

template<typename T>
ompBLAS_status ger_impl(ompBLAS_handle& handle,
                        const int&      m,
                        const int&      n,
                        const T&        alpha,
                        const T* const  x,
                        const int&      incx,
                        const T* const  y,
                        const int&      incy,
                        T* const        A,
                        const int&      lda)
{
  BLAS::ger(m, n, alpha, x, incx, y, incy, A, lda);
  return 0;
}


ompBLAS_status ger(ompBLAS_handle&    handle,
                   const int&         m,
                   const int&         n,
                   const float&       alpha,
                   const float* const x,
                   const int&         incx,
                   const float* const y,
                   const int&         incy,
                   float* const       A,
                   const int&         lda)
{
  return ger_impl(handle, m, n, alpha, x, incx, y, incy, A, lda);
}

ompBLAS_status ger(ompBLAS_handle&     handle,
                   const int&          m,
                   const int&          n,
                   const double&       alpha,
                   const double* const x,
                   const int&          incx,
                   const double* const y,
                   const int&          incy,
                   double* const       A,
                   const int&          lda)
{
  return ger_impl(handle, m, n, alpha, x, incx, y, incy, A, lda);
}

ompBLAS_status ger(ompBLAS_handle&                  handle,
                   const int&                       m,
                   const int&                       n,
                   const std::complex<float>&       alpha,
                   const std::complex<float>* const x,
                   const int&                       incx,
                   const std::complex<float>* const y,
                   const int&                       incy,
                   std::complex<float>* const       A,
                   const int&                       lda)
{
  return ger_impl(handle, m, n, alpha, x, incx, y, incy, A, lda);
}

ompBLAS_status ger(ompBLAS_handle&                   handle,
                   const int&                        m,
                   const int&                        n,
                   const std::complex<double>&       alpha,
                   const std::complex<double>* const x,
                   const int&                        incx,
                   const std::complex<double>* const y,
                   const int&                        incy,
                   std::complex<double>* const       A,
                   const int&                        lda)
{
  return ger_impl(handle, m, n, alpha, x, incx, y, incy, A, lda);
}

} // namespace ompBLAS
} // namespace qmcplusplus
