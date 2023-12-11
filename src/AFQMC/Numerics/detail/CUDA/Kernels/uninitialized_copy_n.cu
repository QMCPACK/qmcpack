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
#include <thrust/complex.h>
#define ENABLE_CUDA 1
#include "AFQMC/Memory/CUDA/cuda_utilities.h"
#include "AFQMC/Numerics/detail/CUDA/Kernels/cuda_settings.h"
//#include "AFQMC/Numerics/detail/CUDA/Kernels/strided_range.hpp"

namespace kernels
{
template<typename T, typename Size>
__global__ void kernel_uninitialized_copy_n(Size N, T const* x, Size incx, T* arr, Size incy)
{
  Size N0(8 * blockDim.x * blockIdx.x);
  T const* x_(x + Size(incx) * N0);
  T* arr_(arr + Size(incy) * N0);
  Size N_(min(Size(8 * blockDim.x), N - N0));
  for (Size ip = Size(threadIdx.x); ip < N_; ip += Size(blockDim.x))
  {
    arr_[ip * incy] = x_[ip * incx];
  }
}


void uninitialized_copy_n(int N, double const* first, int incx, double* array, int incy)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  kernel_uninitialized_copy_n<<<nblk, nthr>>>(N, first, incx, array, incy);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_copy_n(int N, std::complex<double> const* first, int incx, std::complex<double>* array, int incy)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  kernel_uninitialized_copy_n<<<nblk, nthr>>>(N, first, incx, array, incy);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_copy_n(int N, int const* first, int incx, int* array, int incy)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  kernel_uninitialized_copy_n<<<nblk, nthr>>>(N, first, incx, array, incy);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

// long
void uninitialized_copy_n(long N, double const* first, long incx, double* array, long incy)
{
  long N_(8 * DEFAULT_BLOCK_SIZE);
  long nblk((N + N_ - 1) / N_);
  long nthr(DEFAULT_BLOCK_SIZE);
  kernel_uninitialized_copy_n<<<nblk, nthr>>>(N, first, incx, array, incy);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_copy_n(long N, std::complex<double> const* first, long incx, std::complex<double>* array, long incy)
{
  long N_(8 * DEFAULT_BLOCK_SIZE);
  long nblk((N + N_ - 1) / N_);
  long nthr(DEFAULT_BLOCK_SIZE);
  kernel_uninitialized_copy_n<<<nblk, nthr>>>(N, first, incx, array, incy);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_copy_n(long N, int const* first, long incx, int* array, long incy)
{
  long N_(8 * DEFAULT_BLOCK_SIZE);
  long nblk((N + N_ - 1) / N_);
  long nthr(DEFAULT_BLOCK_SIZE);
  kernel_uninitialized_copy_n<<<nblk, nthr>>>(N, first, incx, array, incy);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

} // namespace kernels
