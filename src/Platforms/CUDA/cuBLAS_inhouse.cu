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


#include "Platforms/CUDA/cuBLAS_inhouse.hpp"
#include <stdexcept>
#include <cuComplex.h>

namespace qmcplusplus
{
namespace cuBLAS_inhouse
{
template<typename T, int ROWBS, int COLBS>
__global__ void gemvT_batched_kernel(const int m,
                                     const int n,
                                     const T* __restrict__ alpha,
                                     const T* const A[],
                                     const int lda,
                                     const T* const x[],
                                     const T* __restrict__ beta,
                                     T* const y[])
{
  static_assert(ROWBS <= COLBS, "Row block size must not be larger than column block size!");

  __shared__ T sum[ROWBS][COLBS];
  __shared__ T x_part[COLBS];

  const int tid = threadIdx.x;
  for (int i = 0; i < ROWBS; i++)
    sum[i][tid] = T(0.0);

  const T* __restrict__ A_iw = A[blockIdx.x];
  const T* __restrict__ x_iw = x[blockIdx.x];

  const int row_begin = blockIdx.y * ROWBS;
  const int row_max   = (m - row_begin) < ROWBS ? (m - row_begin) : ROWBS;

  const int num_col_blocks = (n + COLBS - 1) / COLBS;
  for (int ib = 0; ib < num_col_blocks; ib++)
  {
    const int col_id = ib * COLBS + tid;
    if (col_id < n)
      x_part[tid] = x_iw[col_id];
    for (int row_id = row_begin; row_id < row_begin + row_max; row_id++)
      if (col_id < n)
        sum[row_id - row_begin][tid] += x_part[tid] * A_iw[row_id * lda + col_id];
  }
  __syncthreads();

  const int col_max = COLBS < n ? COLBS : n;

  T dot_sum(0);
  for (int col_id = 0; col_id < col_max; col_id++)
    if (tid < row_max)
      dot_sum += sum[tid][col_id];

  T* __restrict__ y_iw = y[blockIdx.x];
  if (tid < row_max)
    y_iw[row_begin + tid] = alpha[blockIdx.x] * dot_sum + beta[blockIdx.x] * y_iw[row_begin + tid];
}

template<typename T>
cuBLAS_inhouse_status gemv_batched_impl(cuBLAS_inhouse_handle& handle,
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
  if (batch_count == 0)
    return cudaSuccess;

  if (trans == 'T')
  {
    if (incx != 1 || incy != 1)
      throw std::runtime_error("incx !=1 or incy != 1 are not implemented in cuBLAS_inhouse::gemv_batched_impl!");

    const int ROWBS          = 32;
    const int COLBS          = 32;
    const int num_row_blocks = (m + ROWBS - 1) / ROWBS;
    dim3 dimBlock(COLBS);
    dim3 dimGrid(batch_count, num_row_blocks);
    gemvT_batched_kernel<T, ROWBS, COLBS><<<dimGrid, dimBlock, 0, handle>>>(m, n, alpha, A, lda, x, beta, y);
    return cudaPeekAtLastError();
  }
  else
  {
    throw std::runtime_error("trans = 'N' not implemented in gemv_batched_impl!");
  }
}

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
                                   const int batch_count)
{
  return gemv_batched_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy, batch_count);
}

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
                                   const int batch_count)
{
  return gemv_batched_impl(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy, batch_count);
}

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
                                   const int batch_count)
{
  //return gemv_batched_impl(handle, trans, m, n, (const cuComplex*)alpha, (const cuComplex* const*)A, lda, (const cuComplex* const*)x, incx, (const cuComplex*)beta, (cuComplex* const*)y, incy, batch_count);
  return cudaSuccess;
}

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
                                   const int batch_count)
{
  //return gemv_batched_impl(handle, trans, m, n, (const cuDoubleComplex*)alpha, (const cuDoubleComplex* const*)A, lda, (const cuDoubleComplex* const*)x, incx, (const cuDoubleComplex*)beta, (cuDoubleComplex* const*)y, incy, batch_count);
  return cudaSuccess;
}


template<typename T>
cuBLAS_inhouse_status ger_batched_impl(cuBLAS_inhouse_handle& handle,
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
  if (batch_count == 0)
    return cudaSuccess;

  if (incx != 1 || incy != 1)
    throw std::runtime_error("incx !=1 or incy != 1 are not implemented in cuBLAS_inhouse::ger_batched_impl!");

  /*
  PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(3) is_device_ptr(A, x, y, alpha)")
  for(size_t ib = 0; ib < batch_count; ib++)
    for(size_t i = 0; i < n; i++)
      for(size_t j = 0; j < m; j++)
        A[ib][i * lda + j] += alpha[ib] * x[ib][j] * y[ib][i];
*/
  return cudaPeekAtLastError();
}

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
                                  const int batch_count)
{
  return ger_batched_impl(handle, m, n, alpha, x, incx, y, incy, A, lda, batch_count);
}

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
                                  const int batch_count)
{
  return ger_batched_impl(handle, m, n, alpha, x, incx, y, incy, A, lda, batch_count);
}

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
                                  const int batch_count)
{
  return ger_batched_impl(handle, m, n, alpha, x, incx, y, incy, A, lda, batch_count);
}

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
                                  const int batch_count)
{
  return ger_batched_impl(handle, m, n, alpha, x, incx, y, incy, A, lda, batch_count);
}
} // namespace cuBLAS_inhouse
} // namespace qmcplusplus
