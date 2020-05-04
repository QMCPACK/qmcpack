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
__global__ void gemvT_batched_kernel(const int m, // number of columns in row major A
                                     const int n, // number of rows in row major A
                                     const T* __restrict__ alpha,
                                     const T* const A[],
                                     const int lda,
                                     const T* const x[],
                                     const int incx,
                                     const T* __restrict__ beta,
                                     T* const y[],
                                     const int incy)
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
  const int row_max   = (n - row_begin) < ROWBS ? (n - row_begin) : ROWBS;

  const int num_col_blocks = (m + COLBS - 1) / COLBS;
  for (int ib = 0; ib < num_col_blocks; ib++)
  {
    const int col_id = ib * COLBS + tid;
    if (col_id < m)
      x_part[tid] = x_iw[col_id * incx];
    for (int row_id = row_begin; row_id < row_begin + row_max; row_id++)
      if (col_id < m)
        sum[row_id - row_begin][tid] += x_part[tid] * A_iw[row_id * lda + col_id];
  }

  for (int iend = COLBS / 2; iend > 0; iend /= 2)
  {
    __syncthreads();
    for (int irow = 0; irow < row_max; irow++)
      if (tid < iend)
        sum[irow][tid] += sum[irow][tid + iend];
  }

  __syncthreads();
  T* __restrict__ y_iw = y[blockIdx.x];
  if (tid < row_max)
    y_iw[(row_begin + tid) * incy] = alpha[blockIdx.x] * sum[tid][0] + beta[blockIdx.x] * y_iw[(row_begin + tid) * incy];
}

template<typename T, int ROWBS>
__global__ void gemvN_batched_kernel(const int m, // number of columns in row major A
                                     const int n, // number of rows in row major A
                                     const T* __restrict__ alpha,
                                     const T* const A[],
                                     const int lda,
                                     const T* const x[],
                                     const int incx,
                                     const T* __restrict__ beta,
                                     T* const y[],
                                     const int incy)
{
  const T* __restrict__ A_iw = A[blockIdx.x];
  const T* __restrict__ x_iw = x[blockIdx.x];
  T* __restrict__ y_iw = y[blockIdx.x];

  const int tid = threadIdx.x;
  const int row_begin = blockIdx.y * ROWBS;

  T sum(0);
  for (int col_id = 0; col_id < n; col_id++)
    sum += A_iw[col_id * lda + row_begin + tid] * x_iw[col_id * incx];
  y_iw[(row_begin + tid) * incy] = alpha[blockIdx.x] * sum + beta[blockIdx.x] * y_iw[(row_begin + tid) * incy];
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
  if (batch_count == 0 || m == 0 || n == 0)
    return cudaSuccess;

  if (trans == 'T')
  {
    const int ROWBS          = 4;
    const int COLBS          = 64;
    const int num_row_blocks = (n + ROWBS - 1) / ROWBS;
    dim3 dimBlock(COLBS);
    dim3 dimGrid(batch_count, num_row_blocks);
    gemvT_batched_kernel<T, ROWBS, COLBS><<<dimGrid, dimBlock, 0, handle>>>(m, n, alpha, A, lda, x, incx, beta, y, incy);
  }
  else
  {
    const int ROWBS          = 64;
    const int num_row_blocks = (m + ROWBS - 1) / ROWBS;
    dim3 dimBlock(ROWBS);
    dim3 dimGrid(batch_count, num_row_blocks);
    gemvN_batched_kernel<T, ROWBS><<<dimGrid, dimBlock, 0, handle>>>(m, n, alpha, A, lda, x, incx, beta, y, incy);
  }
  return cudaPeekAtLastError();
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


template<typename T, int ROWBS, int COLBS>
__global__ void ger_batched_kernel(const int m, // number of columns in row major A
                                   const int n, // number of rows in row major A
                                   const T* __restrict__ alpha,
                                   const T* const x[],
                                   const int incx,
                                   const T* const y[],
                                   const int incy,
                                   T* const A[],
                                   const int lda)
{
  const int iw = blockIdx.x;
  const T* __restrict__ x_iw = x[iw];
  const T* __restrict__ y_iw = y[iw];
  T* __restrict__ A_iw = A[iw];

  const int row_begin = blockIdx.y * ROWBS;
  const int row_end   = (row_begin + ROWBS) < n ? (row_begin + ROWBS) : n;
  const int tid = threadIdx.x;
  const int col_id = blockIdx.z * COLBS + tid;

  __shared__ T x_part[COLBS];
  if (col_id < m)
    x_part[tid] = x_iw[col_id * incx];

  for (int row_id = row_begin; row_id < row_end; row_id++)
    if (col_id < m)
      A_iw[row_id * lda + col_id] += alpha[iw] * x_part[tid] * y_iw[row_id * incy];
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
  if (batch_count == 0 || m == 0 || n == 0)
    return cudaSuccess;

  const int ROWBS          = 16;
  const int COLBS          = 64;
  const int num_row_blocks = (n + ROWBS - 1) / ROWBS;
  const int num_col_blocks = (m + COLBS - 1) / COLBS;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_count, num_row_blocks, num_col_blocks);
  ger_batched_kernel<T, ROWBS, COLBS><<<dimGrid, dimBlock, 0, handle>>>(m, n, alpha, x, incx, y, incy, A, lda);
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
  //return ger_batched_impl(handle, m, n, alpha, x, incx, y, incy, A, lda, batch_count);
  return cudaSuccess;
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
  //return ger_batched_impl(handle, m, n, alpha, x, incx, y, incy, A, lda, batch_count);
  return cudaSuccess;
}

template<typename T, int COLBS>
__global__ void copy_batched_kernel(const int n,
                                    const T* const in[],
                                    T* const out[])
{
  const int iw                 = blockIdx.x;
  const T* __restrict__ in_iw  = in[iw];
  T* __restrict__ out_iw       = out[iw];

  const int col_id = blockIdx.y * COLBS + threadIdx.x;
  out_iw[col_id] = in_iw[col_id];
}

template<typename T>
cuBLAS_inhouse_status copy_batched_impl(cudaStream_t& hstream,
                                   const int n,
                                   const T* const in[],
                                   const int incx,
                                   T* const out[],
                                   const int incy,
                                   const int batch_count)
{
  if (batch_count == 0 || n == 0)
    return cudaSuccess;

  if (incx != 1 || incy != 1)
    throw std::runtime_error("incx !=1 or incy != 1 are not implemented in cuBLAS_inhouse::copy_batched_impl!");

  const int COLBS = 128;
  const int num_col_blocks = (n + COLBS - 1) / COLBS;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_count, num_col_blocks);
  copy_batched_kernel<T, COLBS><<<dimGrid, dimBlock, 0, hstream>>>(n, in, out);
  return cudaPeekAtLastError();
}

cuBLAS_inhouse_status copy_batched(cudaStream_t& hstream,
                                   const int n,
                                   const float* const in[],
                                   const int incx,
                                   float* const out[],
                                   const int incy,
                                   const int batch_count)
{
  return copy_batched_impl(hstream, n, in, incx, out, incy, batch_count);
}

cuBLAS_inhouse_status copy_batched(cudaStream_t& hstream,
                                   const int n,
                                   const double* const in[],
                                   const int incx,
                                   double* const out[],
                                   const int incy,
                                   const int batch_count)
{
  return copy_batched_impl(hstream, n, in, incx, out, incy, batch_count);
}

cuBLAS_inhouse_status copy_batched(cudaStream_t& hstream,
                                   const int n,
                                   const std::complex<float>* const in[],
                                   const int incx,
                                   std::complex<float>* const out[],
                                   const int incy,
                                   const int batch_count)
{
  //return copy_batched_impl(hstream, n, in, incx, out, incy, batch_count);
  return cudaSuccess;
}

cuBLAS_inhouse_status copy_batched(cudaStream_t& hstream,
                                   const int n,
                                   const std::complex<double>* const in[],
                                   const int incx,
                                   std::complex<double>* const out[],
                                   const int incy,
                                   const int batch_count)
{
  //return copy_batched_impl(hstream, n, in, incx, out, incy, batch_count);
  return cudaSuccess;
}

// Xscal_batched
template<typename T, int COLBS>
__global__ void scal_batched_kernel(const int n,
                                    const T* const scale,
                                    T* const inout[],
                                    const int inc)
{
  const int iw             = blockIdx.x;
  T* __restrict__ inout_iw = inout[iw];
  const T scale_iw = scale[iw];

  const int col_id = blockIdx.y * COLBS + threadIdx.x;
  inout_iw[col_id * inc] *= scale_iw;
}

template<typename T>
cuBLAS_inhouse_status scal_batched_impl(cudaStream_t& hstream,
                                   const int n,
                                   const T* const scale,
                                   T* const inout[],
                                   const int inc,
                                   const int batch_count)
{
  if (batch_count == 0 || n == 0)
    return cudaSuccess;

  const int COLBS = 128;
  const int num_col_blocks = (n + COLBS - 1) / COLBS;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_count, num_col_blocks);
  scal_batched_kernel<T, COLBS><<<dimGrid, dimBlock, 0, hstream>>>(n, scale, inout, inc);
  return cudaPeekAtLastError();
}

cuBLAS_inhouse_status scal_batched(cudaStream_t& hstream,
                                   const int n,
                                   const float* const scale,
                                   float* const inout[],
                                   const int inc,
                                   const int batch_count)
{
  return scal_batched_impl(hstream, n, scale, inout, inc, batch_count);
}

cuBLAS_inhouse_status scal_batched(cudaStream_t& hstream,
                                   const int n,
                                   const double* const scale,
                                   double* const inout[],
                                   const int inc,
                                   const int batch_count)
{
  return scal_batched_impl(hstream, n, scale, inout, inc, batch_count);
}

cuBLAS_inhouse_status scal_batched(cudaStream_t& hstream,
                                   const int n,
                                   const std::complex<float>* const scale,
                                   std::complex<float>* const inout[],
                                   const int inc,
                                   const int batch_count)
{
  //return scal_batched_impl(hstream, n, scale, inout, inc, batch_count);
  return cudaSuccess;
}

cuBLAS_inhouse_status scal_batched(cudaStream_t& hstream,
                                   const int n,
                                   const std::complex<double>* const scale,
                                   std::complex<double>* const inout[],
                                   const int inc,
                                   const int batch_count)
{
  //return scal_batched_impl(hstream, n, scale, inout, inc, batch_count);
  return cudaSuccess;
}

// Xdot_batched
template<typename T, int COLBS>
__global__ void dot_batched_kernel(const int n,
                                   const T* const x[],
                                   const int incx,
                                   const T* const y[],
                                   const int incy,
                                   T* const res)
{
  const int tid            = threadIdx.x;
  const int iw             = blockIdx.x;
  const T* __restrict__ x_iw = x[iw];
  const T* __restrict__ y_iw = y[iw];

  __shared__ T sum[COLBS];
  sum[tid] = T(0);
  const int num_col_blocks = (n + COLBS - 1) / COLBS;
  for (int ib = 0; ib < num_col_blocks; ib++)
  {
    const int col_id = ib * COLBS + tid;
    sum[col_id] += x_iw[col_id * incx] * y_iw[col_id * incy];
  }

  for (int iend = COLBS / 2; iend > 0; iend /= 2)
  {
    __syncthreads();
    if (tid < iend)
      sum[tid] += sum[tid + iend];
  }

  if (tid == 0)
    res[iw] = sum[0];
}

template<typename T>
cuBLAS_inhouse_status dot_batched_impl(cudaStream_t& hstream,
                                       const int n,
                                       const T* const x[],
                                       const int incx,
                                       const T* const y[],
                                       const int incy,
                                       T* const res,
                                       const int batch_count)
{
  if (batch_count == 0 || n == 0)
    return cudaSuccess;

  const int COLBS = 128;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_count);
  dot_batched_kernel<T, COLBS><<<dimGrid, dimBlock, 0, hstream>>>(n, x, incx, y, incy, res);
  return cudaPeekAtLastError();
}

cuBLAS_inhouse_status dot_batched(cudaStream_t& hstream,
                                  const int n,
                                  const float* const x[],
                                  const int incx,
                                  const float* const y[],
                                  const int incy,
                                  float* const res,
                                  const int batch_count)
{
  return dot_batched_impl(hstream, n, x, incx, y, incy, res, batch_count);
}

cuBLAS_inhouse_status dot_batched(cudaStream_t& hstream,
                                  const int n,
                                  const double* const x[],
                                  const int incx,
                                  const double* const y[],
                                  const int incy,
                                  double* const res,
                                  const int batch_count)
{
  return dot_batched_impl(hstream, n, x, incx, y, incy, res, batch_count);
}

cuBLAS_inhouse_status dot_batched(cudaStream_t& hstream,
                                  const int n,
                                  const std::complex<float>* const x[],
                                  const int incx,
                                  const std::complex<float>* const y[],
                                  const int incy,
                                  std::complex<float>* const res,
                                  const int batch_count)
{
  //return dot_batched_impl(hstream, n, x, incx, y, incy, res, batch_count);
  return cudaSuccess;
}

cuBLAS_inhouse_status dot_batched(cudaStream_t& hstream,
                                  const int n,
                                  const std::complex<double>* const x[],
                                  const int incx,
                                  const std::complex<double>* const y[],
                                  const int incy,
                                  std::complex<double>* const res,
                                  const int batch_count)
{
  //return dot_batched_impl(hstream, n, x, incx, y, incy, res, batch_count);
  return cudaSuccess;
}

} // namespace cuBLAS_inhouse
} // namespace qmcplusplus
