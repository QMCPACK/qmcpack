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

#define BOOST_NO_AUTO_PTR
#include <cassert>
#include <complex>
#include <type_traits>
#define ENABLE_CUDA 1
#include "AFQMC/Numerics/detail/CUDA/Kernels/cuda_settings.h"
#include "AFQMC/Memory/CUDA/cuda_utilities.h"

namespace kernels
{

template<typename T>
__global__ void kernel_fill_n(long N, T* x, long incx, T const a)
{
  long N0(8 * blockDim.x * blockIdx.x);
  T* x_(x + incx * N0);
  long N_(min(long(8 * blockDim.x), N - N0));
  for (long ip = threadIdx.x; ip < N_; ip += blockDim.x)
  {
    x_[ip * incx] = a;
  }
}

template<typename T>
__global__ void kernel_fill2D_n(long N, long M, T* y, long lda, T const a)
{
  for (long ip = threadIdx.x; ip < N; ip += blockDim.x)
    for (long jp = threadIdx.y; jp < M; jp += blockDim.y)
    {
      y[ip * lda + jp] = a;
    }
}

template<typename T>
void fill_n(T* first, long N, long incx, T const value)
{
  long N_(8 * DEFAULT_BLOCK_SIZE);
  long nblk((N + N_ - 1) / N_);
  long nthr(DEFAULT_BLOCK_SIZE);
  kernel_fill_n<<<nblk, nthr>>>(N, first, incx, value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

template<typename T>
void fill_n(T* first, long N, T const value)
{
  long N_(8 * DEFAULT_BLOCK_SIZE);
  long nblk((N + N_ - 1) / N_);
  long nthr(DEFAULT_BLOCK_SIZE);
  kernel_fill_n<<<nblk, nthr>>>(N, first, 1l, value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

template<typename T>
void fill2D_n(long N, long M, T* A, long lda, T const value)
{
  kernel_fill2D_n<<<32, 32>>>(N, M, A, lda, value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

// template specializations
template void fill_n(char* first, long N, long stride, char const value);
template void fill_n(int* first, long N, long stride, int const value);
template void fill_n(unsigned int* first, long N, long stride, unsigned int const value);
template void fill_n(long* first, long N, long stride, long const value);
template void fill_n(unsigned long* first, long N, long stride, unsigned long const value);
template void fill_n(float* first, long N, long stride, float const value);
template void fill_n(double* first, long N, long stride, double const value);
template void fill_n(std::complex<float>* first, long N, long stride, std::complex<float> const value);
template void fill_n(std::complex<double>* first, long N, long stride, std::complex<double> const value);

template void fill_n(char* first, long N, char const value);
template void fill_n(int* first, long N, int const value);
template void fill_n(unsigned int* first, long N, unsigned int const value);
template void fill_n(long* first, long N, long const value);
template void fill_n(unsigned long* first, long N, unsigned long const value);
template void fill_n(float* first, long N, float const value);
template void fill_n(double* first, long N, double const value);
template void fill_n(std::complex<float>* first, long N, std::complex<float> const value);
template void fill_n(std::complex<double>* first, long N, std::complex<double> const value);

template void fill2D_n(long N, long M, int* A, long lda, int const value);
template void fill2D_n(long N, long M, long* A, long lda, long const value);
template void fill2D_n(long N, long M, unsigned int* A, long lda, unsigned int const value);
template void fill2D_n(long N, long M, unsigned long* A, long lda, unsigned long const value);
template void fill2D_n(long N, long M, float* A, long lda, float const value);
template void fill2D_n(long N, long M, double* A, long lda, double const value);
template void fill2D_n(long N, long M, std::complex<double>* A, long lda, std::complex<double> const value);
template void fill2D_n(long N, long M, std::complex<float>* A, long lda, std::complex<float> const value);

} // namespace kernels

