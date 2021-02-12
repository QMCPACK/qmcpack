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
#include <type_traits>
/*
#include <thrust/complex.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/uninitialized_fill.h>
*/
#define ENABLE_CUDA 1
#include "AFQMC/Memory/CUDA/cuda_utilities.h"
#include "AFQMC/Numerics/detail/CUDA/Kernels/cuda_settings.h"

namespace kernels
{
template<typename T, typename Size>
__global__ void kernel_uninitialized_fill_n(Size N, T* x, T const a)
{
  Size N0(8 * blockDim.x * blockIdx.x);
  T* x_(x + N0);
  Size N_(min(Size(8 * blockDim.x), N - N0));
  for (Size ip = Size(threadIdx.x); ip < N_; ip += Size(blockDim.x))
  {
    x_[ip] = a;
  }
}

void uninitialized_fill_n(bool* first, int N, bool const value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  kernel_uninitialized_fill_n<<<nblk, nthr>>>(N, first, value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_fill_n(int* first, int N, int const value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  kernel_uninitialized_fill_n<<<nblk, nthr>>>(N, first, value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_fill_n(float* first, int N, float const value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  kernel_uninitialized_fill_n<<<nblk, nthr>>>(N, first, value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_fill_n(double* first, int N, double const value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  kernel_uninitialized_fill_n<<<nblk, nthr>>>(N, first, value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_fill_n(std::complex<float>* first, int N, std::complex<float> const value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  kernel_uninitialized_fill_n<<<nblk, nthr>>>(N, first, value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_fill_n(std::complex<double>* first, int N, std::complex<double> const value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  kernel_uninitialized_fill_n<<<nblk, nthr>>>(N, first, value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

/*void uninitialized_fill_n(double2 * first, int N, double2 const value)*/
/*{ */
/*kernel_uninitialized_fill_n<<<1,256>>>(N,first,value);*/
/*qmc_cuda::cuda_check(cudaGetLastError());*/
/*qmc_cuda::cuda_check(cudaDeviceSynchronize());*/
/*}*/

void uninitialized_fill_n(bool* first, long N, bool const value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  kernel_uninitialized_fill_n<<<nblk, nthr>>>(N, first, value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_fill_n(int* first, long N, int const value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  kernel_uninitialized_fill_n<<<nblk, nthr>>>(N, first, value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_fill_n(float* first, long N, float const value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  kernel_uninitialized_fill_n<<<nblk, nthr>>>(N, first, value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_fill_n(double* first, long N, double const value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  kernel_uninitialized_fill_n<<<nblk, nthr>>>(N, first, value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_fill_n(std::complex<float>* first, long N, std::complex<float> const value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  kernel_uninitialized_fill_n<<<nblk, nthr>>>(N, first, value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_fill_n(std::complex<double>* first, long N, std::complex<double> const value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  kernel_uninitialized_fill_n<<<nblk, nthr>>>(N, first, value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

/*void uninitialized_fill_n(double2 * first, long N, double2 const value)*/
/*{*/
/*kernel_uninitialized_fill_n<<<1,256>>>(N,first,value);*/
/*qmc_cuda::cuda_check(cudaGetLastError());*/
/*qmc_cuda::cuda_check(cudaDeviceSynchronize());*/
/*}*/


} // namespace kernels
