//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//    Lawrence Livermore National Laboratory
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <complex>
#include <cuda.h>
#include <thrust/complex.h>
#include <cuda_runtime.h>
#include "AFQMC/Numerics/detail/CUDA/Kernels/cuda_settings.h"
#define ENABLE_CUDA 1
#include "AFQMC/Memory/CUDA/cuda_utilities.h"

namespace kernels
{
// many ways to do this! sloppy for now!
template<typename T, typename T2>
__global__ void kernel_tbt_mv_plus(int dim, int nrow, int ncol, T* A, int lda, T2 const* x, int incx)
{
  // assuming a single block
  if (blockIdx.x > 0 || blockIdx.y > 0 || blockIdx.z > 0)
    return;
  if (threadIdx.x > nrow)
    return;
  if (threadIdx.y > ncol)
    return;
  if (dim == 0)
  {
    for (int i = threadIdx.x; i < nrow; i += blockDim.x)
      for (int j = threadIdx.y; j < ncol; j += blockDim.y)
        A[i * lda + j] += x[i * incx];
  }
  else
  {
    for (int i = threadIdx.x; i < nrow; i += blockDim.x)
      for (int j = threadIdx.y; j < ncol; j += blockDim.y)
        A[i * lda + j] += x[j * incx];
  }
}

template<typename T, typename T2>
__global__ void kernel_tbt_mv_minus(int dim, int nrow, int ncol, T* A, int lda, T2 const* x, int incx)
{
  // assuming a single block
  if (blockIdx.x > 0 || blockIdx.y > 0 || blockIdx.z > 0)
    return;
  if (threadIdx.x > nrow)
    return;
  if (threadIdx.y > ncol)
    return;
  if (dim == 0)
  {
    for (int i = threadIdx.x; i < nrow; i += blockDim.x)
      for (int j = threadIdx.y; j < ncol; j += blockDim.y)
        A[i * lda + j] -= x[i * incx];
  }
  else
  {
    for (int i = threadIdx.x; i < nrow; i += blockDim.x)
      for (int j = threadIdx.y; j < ncol; j += blockDim.y)
        A[i * lda + j] -= x[j * incx];
  }
}

template<typename T, typename T2>
__global__ void kernel_tbt_mv_mult(int dim, int nrow, int ncol, T* A, int lda, T2 const* x, int incx)
{
  // assuming a single block
  if (blockIdx.x > 0 || blockIdx.y > 0 || blockIdx.z > 0)
    return;
  if (threadIdx.x > nrow)
    return;
  if (threadIdx.y > ncol)
    return;
  if (dim == 0)
  {
    for (int i = threadIdx.x; i < nrow; i += blockDim.x)
      for (int j = threadIdx.y; j < ncol; j += blockDim.y)
        A[i * lda + j] *= x[i * incx];
  }
  else
  {
    for (int i = threadIdx.x; i < nrow; i += blockDim.x)
      for (int j = threadIdx.y; j < ncol; j += blockDim.y)
        A[i * lda + j] *= x[j * incx];
  }
}

template<typename T, typename T2>
__global__ void kernel_tbt_mv_div(int dim, int nrow, int ncol, T* A, int lda, T2 const* x, int incx)
{
  // assuming a single block
  if (blockIdx.x > 0 || blockIdx.y > 0 || blockIdx.z > 0)
    return;
  if (threadIdx.x > nrow)
    return;
  if (threadIdx.y > ncol)
    return;
  if (dim == 0)
  {
    for (int i = threadIdx.x; i < nrow; i += blockDim.x)
      for (int j = threadIdx.y; j < ncol; j += blockDim.y)
        A[i * lda + j] /= x[i * incx];
  }
  else
  {
    for (int i = threadIdx.x; i < nrow; i += blockDim.x)
      for (int j = threadIdx.y; j < ncol; j += blockDim.y)
        A[i * lda + j] /= x[j * incx];
  }
}

void term_by_term_mat_vec_plus(int dim,
                               int nrow,
                               int ncol,
                               std::complex<double>* A,
                               int lda,
                               std::complex<double>* x,
                               int incx)
{
  int xblock_dim = std::min(nrow, 32);
  int yblock_dim = std::min(ncol, 32);
  dim3 block_dim(xblock_dim, yblock_dim, 1);
  kernel_tbt_mv_plus<<<1, block_dim>>>(dim, nrow, ncol, reinterpret_cast<thrust::complex<double>*>(A), lda,
                                       reinterpret_cast<thrust::complex<double> const*>(x), incx);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void term_by_term_mat_vec_minus(int dim,
                                int nrow,
                                int ncol,
                                std::complex<double>* A,
                                int lda,
                                std::complex<double>* x,
                                int incx)
{
  int xblock_dim = std::min(nrow, 32);
  int yblock_dim = std::min(ncol, 32);
  dim3 block_dim(xblock_dim, yblock_dim, 1);
  kernel_tbt_mv_minus<<<1, block_dim>>>(dim, nrow, ncol, reinterpret_cast<thrust::complex<double>*>(A), lda,
                                        reinterpret_cast<thrust::complex<double> const*>(x), incx);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void term_by_term_mat_vec_mult(int dim,
                               int nrow,
                               int ncol,
                               std::complex<double>* A,
                               int lda,
                               std::complex<double>* x,
                               int incx)
{
  int xblock_dim = std::min(nrow, 32);
  int yblock_dim = std::min(ncol, 32);
  dim3 block_dim(xblock_dim, yblock_dim, 1);
  kernel_tbt_mv_mult<<<1, block_dim>>>(dim, nrow, ncol, reinterpret_cast<thrust::complex<double>*>(A), lda,
                                       reinterpret_cast<thrust::complex<double> const*>(x), incx);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void term_by_term_mat_vec_div(int dim,
                              int nrow,
                              int ncol,
                              std::complex<double>* A,
                              int lda,
                              std::complex<double>* x,
                              int incx)
{
  int xblock_dim = std::min(nrow, 32);
  int yblock_dim = std::min(ncol, 32);
  dim3 block_dim(xblock_dim, yblock_dim, 1);
  kernel_tbt_mv_div<<<1, block_dim>>>(dim, nrow, ncol, reinterpret_cast<thrust::complex<double>*>(A), lda,
                                      reinterpret_cast<thrust::complex<double> const*>(x), incx);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void term_by_term_mat_vec_plus(int dim, int nrow, int ncol, std::complex<double>* A, int lda, double* x, int incx)
{
  int xblock_dim = std::min(nrow, 32);
  int yblock_dim = std::min(ncol, 32);
  dim3 block_dim(xblock_dim, yblock_dim, 1);
  kernel_tbt_mv_plus<<<1, block_dim>>>(dim, nrow, ncol, reinterpret_cast<thrust::complex<double>*>(A), lda, x, incx);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void term_by_term_mat_vec_minus(int dim, int nrow, int ncol, std::complex<double>* A, int lda, double* x, int incx)
{
  int xblock_dim = std::min(nrow, 32);
  int yblock_dim = std::min(ncol, 32);
  dim3 block_dim(xblock_dim, yblock_dim, 1);
  kernel_tbt_mv_minus<<<1, block_dim>>>(dim, nrow, ncol, reinterpret_cast<thrust::complex<double>*>(A), lda, x, incx);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void term_by_term_mat_vec_mult(int dim, int nrow, int ncol, std::complex<double>* A, int lda, double* x, int incx)
{
  int xblock_dim = std::min(nrow, 32);
  int yblock_dim = std::min(ncol, 32);
  dim3 block_dim(xblock_dim, yblock_dim, 1);
  kernel_tbt_mv_mult<<<1, block_dim>>>(dim, nrow, ncol, reinterpret_cast<thrust::complex<double>*>(A), lda, x, incx);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void term_by_term_mat_vec_div(int dim, int nrow, int ncol, std::complex<double>* A, int lda, double* x, int incx)
{
  int xblock_dim = std::min(nrow, 32);
  int yblock_dim = std::min(ncol, 32);
  dim3 block_dim(xblock_dim, yblock_dim, 1);
  kernel_tbt_mv_div<<<1, block_dim>>>(dim, nrow, ncol, reinterpret_cast<thrust::complex<double>*>(A), lda, x, incx);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}


} // namespace kernels
