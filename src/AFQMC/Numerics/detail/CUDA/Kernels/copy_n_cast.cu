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
#define ENABLE_CUDA 1
#include "AFQMC/Memory/CUDA/cuda_utilities.h"
#include "AFQMC/Numerics/detail/CUDA/Kernels/cuda_settings.h"

namespace kernels
{
template<typename T, typename Q>
__global__ void kernel_copy_n_cast(T const* A, long N, Q* B)
{
  long N0(8 * blockDim.x * blockIdx.x);
  T const* A_(A + N0);
  Q* B_(B + N0);
  long N_(min(long(8 * blockDim.x), N - N0));
  for (long ip = threadIdx.x; ip < N_; ip += blockDim.x)
  {
    B_[ip] = static_cast<Q>(A_[ip]);
  }
}

template<typename T, typename Q>
__global__ void kernel_copy_n_cast(thrust::complex<T> const* A, long N, thrust::complex<Q>* B)
{
  long N0(8 * blockDim.x * blockIdx.x);
  thrust::complex<T> const* A_(A + N0);
  thrust::complex<Q>* B_(B + N0);
  long N_(min(long(8 * blockDim.x), N - N0));
  for (long ip = threadIdx.x; ip < N_; ip += blockDim.x)
  {
    B_[ip] = static_cast<thrust::complex<Q>>(A_[ip]);
  }
}

void copy_n_cast(double const* A, long n, float* B)
{
  long n_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((n + n_ - 1) / n_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  kernel_copy_n_cast<<<nblk, nthr>>>(A, n, B);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void copy_n_cast(float const* A, long n, double* B)
{
  long n_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((n + n_ - 1) / n_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  kernel_copy_n_cast<<<nblk, nthr>>>(A, n, B);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void copy_n_cast(std::complex<double> const* A, long n, std::complex<float>* B)
{
  long n_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((n + n_ - 1) / n_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  kernel_copy_n_cast<<<nblk, nthr>>>(reinterpret_cast<thrust::complex<double> const*>(A), n,
                                     reinterpret_cast<thrust::complex<float>*>(B));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void copy_n_cast(std::complex<float> const* A, long n, std::complex<double>* B)
{
  long n_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((n + n_ - 1) / n_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  kernel_copy_n_cast<<<nblk, nthr>>>(reinterpret_cast<thrust::complex<float> const*>(A), n,
                                     reinterpret_cast<thrust::complex<double>*>(B));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

} // namespace kernels
