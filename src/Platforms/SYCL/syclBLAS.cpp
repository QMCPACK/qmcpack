//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "syclBLAS.hpp"
#include "oneapi/mkl/blas.hpp"

namespace qmcplusplus
{
namespace syclBLAS
{
template<typename T>
sycl::event gemv(sycl::queue& handle,
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
                 const int incy,
                 const std::vector<sycl::event>& events)
{
  return oneapi::mkl::blas::gemv(handle, convertTransEnum(trans), m, n, alpha, A, lda, x, incx, beta, y, incy, events);
}

template sycl::event gemv(sycl::queue& handle,
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
                          const int incy,
                          const std::vector<sycl::event>& events);

template sycl::event gemv(sycl::queue& handle,
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
                          const int incy,
                          const std::vector<sycl::event>& events);

template sycl::event gemv(sycl::queue& handle,
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
                          const int incy,
                          const std::vector<sycl::event>& events);

template sycl::event gemv(sycl::queue& handle,
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
                          const int incy,
                          const std::vector<sycl::event>& events);

/** gemv trans = 'T' case. COLS refers to columns of the m x n column-major Fortran matrix A.
 */
template<typename T, unsigned COLBS>
sycl::event gemvT_batched_impl(sycl::queue& handle,
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
                               const size_t batch_count,
                               const std::vector<sycl::event>& events = {})
{
  if (m == 0 || n == 0 || batch_count == 0)
    return sycl::event();

  // TODO: templatize?
  static constexpr int ROWBS = 4;

  constexpr int SUM_SIZE = ROWBS * COLBS;
  const int num_row_blocks = (n + ROWBS - 1) / ROWBS;

  const size_t num_groups = static_cast<size_t>(batch_count) * static_cast<size_t>(num_row_blocks);
  const sycl::nd_range<1> nd{num_groups * COLBS, COLBS};

  return handle.submit([&](sycl::handler& h) {
    if (!events.empty())
      h.depends_on(events);

    sycl::local_accessor<T, 1> sum(SUM_SIZE, h);

    h.parallel_for(nd, [=](sycl::nd_item<1> item) {
      const int tid     = item.get_local_id(0); // threadIdx.x
      const size_t gid  = item.get_group(0);
      const int by      = static_cast<int>(gid % num_row_blocks); // blockIdx.y
      const int bx      = static_cast<int>(gid / num_row_blocks); // blockIdx.x

      for (int i = 0; i < ROWBS; ++i)
        sum[i * COLBS + tid] = T(0);

      const T* __restrict__ A_iw = A[bx];
      const T* __restrict__ x_iw = x[bx];

      const int row_begin = by * ROWBS;
      const int row_max   = (n - row_begin) < ROWBS ? (n - row_begin) : ROWBS;

      const int num_col_blocks = (m + COLBS - 1) / COLBS;
      for (int ib = 0; ib < num_col_blocks; ++ib)
      {
        if (const int col_id = ib * COLBS + tid; col_id < m)
        {
          const T xv = x_iw[col_id * incx];
          for (int row_id = row_begin; row_id < row_begin + row_max; ++row_id)
            sum[(row_id - row_begin) * COLBS + tid] += xv * A_iw[row_id * lda + col_id];
        }
      }

      for (int iend = COLBS / 2; iend > 0; iend /= 2)
      {
        item.barrier(sycl::access::fence_space::local_space);
        for (int irow = 0; irow < row_max; ++irow)
          if (tid < iend)
            sum[irow * COLBS + tid] += sum[irow * COLBS + tid + iend];
      }

      item.barrier(sycl::access::fence_space::local_space);
      T* __restrict__ y_iw = y[bx];
      if (tid < row_max)
      {
        const int yidx = (row_begin + tid) * incy;
        const T s      = sum[tid * COLBS];
        if (beta[bx] == T(0))
          y_iw[yidx] = alpha[bx] * s; // protect NaN from y_iw
        else
          y_iw[yidx] = alpha[bx] * s + beta[bx] * y_iw[yidx];
      }
    });
  });
}

/** gemv trans = 'N' case. ROW refers to rows of the m x n column-major Fortran matrix A.
 */
template<typename T, unsigned ROWBS>
sycl::event gemvN_batched_impl(sycl::queue& handle,
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
                               const size_t batch_count,
                               const std::vector<sycl::event>& events = {})
{
  if (m == 0 || n == 0 || batch_count == 0)
    return sycl::event();

  const int num_row_blocks = (m + ROWBS - 1) / ROWBS;
  return handle.parallel_for(sycl::nd_range<2>{{batch_count, num_row_blocks * ROWBS}, {1, ROWBS}},
                             [=](sycl::nd_item<2> item) {
                               const unsigned batch = item.get_group(0);
                               const int row        = item.get_global_id(1);
                               if (row < m)
                               {
                                 T sum(0);
                                 for (int col = 0; col < n; col++)
                                   sum += A[batch][col * lda + row] * x[batch][col * incx];
                                 if (beta[batch] == T(0))
                                   y[batch][row * incy] = alpha[batch] * sum; // protecting NaN from y_iw
                                 else
                                   y[batch][row * incy] = alpha[batch] * sum + beta[batch] * y[batch][row * incy];
                               }
                             });
}

template<>
sycl::event gemv_batched<float>(sycl::queue& handle,
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
                                const size_t batch_count,
                                const std::vector<sycl::event>& events)
{
  if (trans == 'N' || trans == 'n')
    return gemvN_batched_impl<float, 64>(handle, m, n, alpha, A, lda, x, incx, beta, y, incy, batch_count);
  else if (trans == 'T' || trans == 't')
    return gemvT_batched_impl<float, 64>(handle, m, n, alpha, A, lda, x, incx, beta, y, incy, batch_count);
  else
    throw std::runtime_error("syclBLAS::gemv_batched only supports 'N', 'T', 'C', 'n'. Input value is " +
                             std::string(1, trans));
}

template<>
sycl::event gemv_batched<double>(sycl::queue& handle,
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
                                 const size_t batch_count,
                                 const std::vector<sycl::event>& events)
{
  if (trans == 'N' || trans == 'n')
    return gemvN_batched_impl<double, 64>(handle, m, n, alpha, A, lda, x, incx, beta, y, incy, batch_count);
  else if (trans == 'T' || trans == 't')
    return gemvT_batched_impl<double, 64>(handle, m, n, alpha, A, lda, x, incx, beta, y, incy, batch_count);
  else
    throw std::runtime_error("syclBLAS::gemv_batched only supports 'N', 'T', 'C', 'n'. Input value is " +
                             std::string(1, trans));
}

template<>
sycl::event gemv_batched<std::complex<float>>(sycl::queue& handle,
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
                                              const size_t batch_count,
                                              const std::vector<sycl::event>& events)
{
  if (trans == 'N' || trans == 'n')
    return gemvN_batched_impl<std::complex<float>, 64>(handle, m, n, alpha, A, lda, x, incx, beta, y, incy,
                                                       batch_count);
  else if (trans == 'T' || trans == 't')
    return gemvT_batched_impl<std::complex<float>, 64>(handle, m, n, alpha, A, lda, x, incx, beta, y, incy,
                                                       batch_count);
  else
    throw std::runtime_error("syclBLAS::gemv_batched only supports 'N', 'T', 'C', 'n'. Input value is " +
                             std::string(1, trans));
}

template<>
sycl::event gemv_batched<std::complex<double>>(sycl::queue& handle,
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
                                               const size_t batch_count,
                                               const std::vector<sycl::event>& events)
{
  if (trans == 'N' || trans == 'n')
    return gemvN_batched_impl<std::complex<double>, 64>(handle, m, n, alpha, A, lda, x, incx, beta, y, incy,
                                                        batch_count);
  else if (trans == 'T' || trans == 't')
    return gemvT_batched_impl<std::complex<double>, 64>(handle, m, n, alpha, A, lda, x, incx, beta, y, incy,
                                                        batch_count);
  else
    throw std::runtime_error("syclBLAS::gemv_batched only supports 'N', 'T', 'C', 'n'. Input value is " +
                             std::string(1, trans));
}

template<typename T>
sycl::event gemm(sycl::queue& handle,
                 const char tA,
                 const char tB,
                 const int m,
                 const int n,
                 const int k,
                 const T alpha,
                 const T* A,
                 const int lda,
                 const T* B,
                 const int ldb,
                 const T beta,
                 T* C,
                 const int ldc,
                 const std::vector<sycl::event>& events)
{
  return oneapi::mkl::blas::gemm(handle, convertTransEnum(tA), convertTransEnum(tB), m, n, k, alpha, A, lda, B, ldb,
                                 beta, C, ldc, events);
}


template sycl::event gemm(sycl::queue& handle,
                          const char tA,
                          const char tB,
                          const int m,
                          const int n,
                          const int k,
                          const float alpha,
                          const float* const A,
                          const int lda,
                          const float* const B,
                          const int ldb,
                          const float beta,
                          float* const C,
                          const int ldc,
                          const std::vector<sycl::event>& events);

template sycl::event gemm(sycl::queue& handle,
                          const char tA,
                          const char tB,
                          const int m,
                          const int n,
                          const int k,
                          const double alpha,
                          const double* const A,
                          const int lda,
                          const double* const B,
                          const int ldb,
                          const double beta,
                          double* const C,
                          const int ldc,
                          const std::vector<sycl::event>& events);

template sycl::event gemm(sycl::queue& handle,
                          const char tA,
                          const char tB,
                          const int m,
                          const int n,
                          const int k,
                          const std::complex<float> alpha,
                          const std::complex<float>* const A,
                          const int lda,
                          const std::complex<float>* const B,
                          const int ldb,
                          const std::complex<float> beta,
                          std::complex<float>* const C,
                          const int ldc,
                          const std::vector<sycl::event>& events);

template sycl::event gemm(sycl::queue& handle,
                          const char tA,
                          const char tB,
                          const int m,
                          const int n,
                          const int k,
                          const std::complex<double> alpha,
                          const std::complex<double>* const A,
                          const int lda,
                          const std::complex<double>* const B,
                          const int ldb,
                          const std::complex<double> beta,
                          std::complex<double>* const C,
                          const int ldc,
                          const std::vector<sycl::event>& events);

template<typename T, int TILE_SIZE, int ROWBS>
sycl::event ger_batched_impl(sycl::queue& handle,
                             const int m,
                             const int n,
                             const T* alpha,
                             const T* const x[],
                             const int incx,
                             const T* const y[],
                             const int incy,
                             T* const A[],
                             const int lda,
                             const size_t batch_count,
                             const std::vector<sycl::event>& events)
{
  static_assert(ROWBS <= TILE_SIZE, "ROWBS cannot be larger than TILE_SIZE!");
  if (m == 0 || n == 0 || batch_count == 0)
    return sycl::event();

  // A is m x n in Fortran, n x m in C.
  constexpr size_t tile_size  = TILE_SIZE;
  constexpr size_t block_rows = ROWBS;
  // the computation is tiled and distributed.
  const size_t row_tiles = (n + tile_size - 1) / tile_size;
  const size_t col_tiles = (m + tile_size - 1) / tile_size;

  return handle.parallel_for(sycl::nd_range<3>{{batch_count, row_tiles * block_rows, col_tiles * tile_size},
                                               {1, block_rows, tile_size}},
                             [=](sycl::nd_item<3> item) {
                               const unsigned batch      = item.get_group(0);
                               const unsigned thX        = item.get_local_id(2);
                               const unsigned thY        = item.get_local_id(1);
                               const unsigned column     = item.get_group(2) * tile_size + thX;
                               const unsigned row_offset = item.get_group(1) * tile_size + thY;
                               if (column < m)
                               {
                                 const T alphaX = alpha[batch] * x[batch][column * incx];
                                 for (unsigned j = 0; j < tile_size; j += block_rows)
                                   if (const unsigned row = row_offset + j; row < n)
                                     A[batch][row * lda + column] += alphaX * y[batch][row * incy];
                               }
                             });
}

template<>
sycl::event ger_batched<float>(sycl::queue& handle,
                               const int m,
                               const int n,
                               const float* alpha,
                               const float* const x[],
                               const int incx,
                               const float* const y[],
                               const int incy,
                               float* const A[],
                               const int lda,
                               const size_t batch_count,
                               const std::vector<sycl::event>& events)
{
  return ger_batched_impl<float, 32, 8>(handle, m, n, alpha, x, incx, y, incy, A, lda, batch_count, events);
}

template<>
sycl::event ger_batched<double>(sycl::queue& handle,
                                const int m,
                                const int n,
                                const double* alpha,
                                const double* const x[],
                                const int incx,
                                const double* const y[],
                                const int incy,
                                double* const A[],
                                const int lda,
                                const size_t batch_count,
                                const std::vector<sycl::event>& events)
{
  return ger_batched_impl<double, 32, 8>(handle, m, n, alpha, x, incx, y, incy, A, lda, batch_count, events);
}

template<>
sycl::event ger_batched<std::complex<float>>(sycl::queue& handle,
                                             const int m,
                                             const int n,
                                             const std::complex<float>* alpha,
                                             const std::complex<float>* const x[],
                                             const int incx,
                                             const std::complex<float>* const y[],
                                             const int incy,
                                             std::complex<float>* const A[],
                                             const int lda,
                                             const size_t batch_count,
                                             const std::vector<sycl::event>& events)
{
  return ger_batched_impl<std::complex<float>, 32, 8>(handle, m, n, alpha, x, incx, y, incy, A, lda, batch_count,
                                                      events);
}

template<>
sycl::event ger_batched<std::complex<double>>(sycl::queue& handle,
                                              const int m,
                                              const int n,
                                              const std::complex<double>* alpha,
                                              const std::complex<double>* const x[],
                                              const int incx,
                                              const std::complex<double>* const y[],
                                              const int incy,
                                              std::complex<double>* const A[],
                                              const int lda,
                                              const size_t batch_count,
                                              const std::vector<sycl::event>& events)
{
  return ger_batched_impl<std::complex<double>, 32, 8>(handle, m, n, alpha, x, incx, y, incy, A, lda, batch_count,
                                                       events);
}

//transpose
template<typename T1, typename T2>
sycl::event transpose(sycl::queue& q,
                      const T1* restrict in,
                      int m,
                      int lda,
                      T2* restrict out,
                      int n,
                      int ldb,
                      const std::vector<sycl::event>& events)
{
  constexpr size_t tile_size = 16;
  const size_t m_max         = ((m + tile_size - 1) / tile_size) * tile_size;
  const size_t n_max         = ((n + tile_size - 1) / tile_size) * tile_size;

  return q.submit([&](sycl::handler& cgh) {
    cgh.depends_on(events);
    sycl::local_accessor<T2, 2> tile(sycl::range<2>(tile_size, tile_size + 1), cgh);

    cgh.parallel_for(sycl::nd_range<2>{{m_max, n_max}, {tile_size, tile_size}}, [=](sycl::nd_item<2> item) {
      unsigned x   = item.get_global_id(1);
      unsigned y   = item.get_global_id(0);
      unsigned xth = item.get_local_id(1);
      unsigned yth = item.get_local_id(0);

      if (x < n && y < m)
        tile[yth][xth] = in[(y)*lda + x];
      item.barrier(sycl::access::fence_space::local_space);

      x = item.get_group(0) * tile_size + xth;
      y = item.get_group(1) * tile_size + yth;
      if (x < m && y < n)
        out[(y)*ldb + x] = tile[xth][yth];
    });
  });
}

template sycl::event transpose(sycl::queue& q,
                               const float* restrict in,
                               int m,
                               int lda,
                               double* restrict out,
                               int n,
                               int ldb,
                               const std::vector<sycl::event>& events);

template sycl::event transpose(sycl::queue& q,
                               const double* restrict in,
                               int m,
                               int lda,
                               double* restrict out,
                               int n,
                               int ldb,
                               const std::vector<sycl::event>& events);

template sycl::event transpose(sycl::queue& q,
                               const std::complex<float>* restrict in,
                               int m,
                               int lda,
                               std::complex<double>* restrict out,
                               int n,
                               int ldb,
                               const std::vector<sycl::event>& events);

template sycl::event transpose(sycl::queue& q,
                               const std::complex<double>* restrict in,
                               int m,
                               int lda,
                               std::complex<double>* restrict out,
                               int n,
                               int ldb,
                               const std::vector<sycl::event>& events);

//copy_n for mixed precision
template<typename T1, typename T2>
sycl::event copy_n(sycl::queue& aq,
                   const T1* restrict VA,
                   size_t array_size,
                   T2* restrict VC,
                   const std::vector<sycl::event>& events)
{
  if (array_size == 0)
    return sycl::event();
  constexpr size_t tile_size = 64;
  const size_t a_max         = ((array_size + tile_size - 1) / tile_size) * tile_size;
  return aq.parallel_for(sycl::range<1>{a_max}, events, [=](sycl::id<1> id) {
    if (id < array_size)
      VC[id] = static_cast<T2>(VA[id]);
  });
}

template sycl::event copy_n(sycl::queue& aq,
                            const double* restrict VA,
                            size_t array_size,
                            float* restrict VC,
                            const std::vector<sycl::event>& events);

template sycl::event copy_n(sycl::queue& aq,
                            const std::complex<double>* restrict VA,
                            size_t array_size,
                            std::complex<float>* restrict VC,
                            const std::vector<sycl::event>& events);

} // namespace syclBLAS

} // namespace qmcplusplus
