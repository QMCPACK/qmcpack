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
#include "AFQMC/Memory/CUDA/cuda_utilities.h"
#include "AFQMC/Numerics/detail/CUDA/Kernels/cuda_settings.h"

namespace kernels
{
template<typename T1, typename T>
__global__ void kernel_inplace_product(int nbatch, int n, int m, T1 const* B, int ldb, T* A, int lda)
{
  int nb = blockIdx.z;
  int i  = blockIdx.x;
  if ((i < n) && (nb < nbatch))
  {
    int j = blockDim.y * blockIdx.y + threadIdx.y;
    if (j < m)
      *(A + (nb * n + i) * lda + j) *= static_cast<T>((*(B + i * ldb + j)));
  }
}

void inplace_product(int nbatch, int n, int m, double const* B, int ldb, std::complex<double>* A, int lda)
{
  int nby = (m + MAX_THREADS_PER_DIM - 1L) / MAX_THREADS_PER_DIM;
  dim3 grid_dim(n, nby, nbatch);
  dim3 block_dim(1, MAX_THREADS_PER_DIM, 1);
  kernel_inplace_product<<<grid_dim, block_dim>>>(nbatch, n, m, B, ldb, reinterpret_cast<thrust::complex<double>*>(A),
                                                  lda);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void inplace_product(int nbatch, int n, int m, std::complex<double> const* B, int ldb, std::complex<double>* A, int lda)
{
  int nby = (m + MAX_THREADS_PER_DIM - 1L) / MAX_THREADS_PER_DIM;
  dim3 grid_dim(n, nby, nbatch);
  dim3 block_dim(1, MAX_THREADS_PER_DIM, 1);
  kernel_inplace_product<<<grid_dim, block_dim>>>(nbatch, n, m, reinterpret_cast<thrust::complex<double> const*>(B),
                                                  ldb, reinterpret_cast<thrust::complex<double>*>(A), lda);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void inplace_product(int nbatch, int n, int m, float const* B, int ldb, std::complex<float>* A, int lda)
{
  int nby = (m + MAX_THREADS_PER_DIM - 1L) / MAX_THREADS_PER_DIM;
  dim3 grid_dim(n, nby, nbatch);
  dim3 block_dim(1, MAX_THREADS_PER_DIM, 1);
  kernel_inplace_product<<<grid_dim, block_dim>>>(nbatch, n, m, B, ldb, reinterpret_cast<thrust::complex<float>*>(A),
                                                  lda);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void inplace_product(int nbatch, int n, int m, std::complex<float> const* B, int ldb, std::complex<float>* A, int lda)
{
  int nby = (m + MAX_THREADS_PER_DIM - 1L) / MAX_THREADS_PER_DIM;
  dim3 grid_dim(n, nby, nbatch);
  dim3 block_dim(1, MAX_THREADS_PER_DIM, 1);
  kernel_inplace_product<<<grid_dim, block_dim>>>(nbatch, n, m, reinterpret_cast<thrust::complex<float> const*>(B), ldb,
                                                  reinterpret_cast<thrust::complex<float>*>(A), lda);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

} // namespace kernels
