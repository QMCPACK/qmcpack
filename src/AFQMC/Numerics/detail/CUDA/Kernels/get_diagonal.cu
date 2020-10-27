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
// simple
// A[k][i] = B[k][i][i]
template<typename T>
__global__ void kernel_get_diagonal_strided(int nk,
                                            int ni,
                                            thrust::complex<T> const* B,
                                            int ldb,
                                            int stride,
                                            thrust::complex<T>* A,
                                            int lda)
{
  int k = blockIdx.y;
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if ((i < ni) && (k < nk))
    A[k * lda + i] = B[k * stride + i * ldb + i];
}

// A[k][i] = B[k][i][i]
void get_diagonal_strided(int nk,
                          int ni,
                          std::complex<double> const* B,
                          int ldb,
                          int stride,
                          std::complex<double>* A,
                          int lda)
{
  size_t nthr = 32;
  size_t nbks = (ni + nthr - 1) / nthr;
  dim3 grid_dim(nbks, nk, 1);
  kernel_get_diagonal_strided<<<grid_dim, nthr>>>(nk, ni, reinterpret_cast<thrust::complex<double> const*>(B), ldb,
                                                  stride, reinterpret_cast<thrust::complex<double>*>(A), lda);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void get_diagonal_strided(int nk,
                          int ni,
                          std::complex<float> const* B,
                          int ldb,
                          int stride,
                          std::complex<float>* A,
                          int lda)
{
  size_t nthr = 32;
  size_t nbks = (ni + nthr - 1) / nthr;
  dim3 grid_dim(nbks, nk, 1);
  kernel_get_diagonal_strided<<<grid_dim, nthr>>>(nk, ni, reinterpret_cast<thrust::complex<float> const*>(B), ldb,
                                                  stride, reinterpret_cast<thrust::complex<float>*>(A), lda);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

} // namespace kernels
