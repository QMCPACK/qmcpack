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
//#include "AFQMC/Numerics/detail/CUDA/Kernels/strided_range.hpp"

namespace kernels
{
template<typename T, typename Size>
__global__ void kernel_uninitialized_copy_n(Size N, T const* x, Size incx, T* arr, Size incy)
{
  for (Size ip = Size(threadIdx.x); ip < N; ip += Size(blockDim.x))
  {
    arr[ip * incy] = x[ip * incx];
  }
}


void uninitialized_copy_n(int N, double const* first, int incx, double* array, int incy)
{
  kernel_uninitialized_copy_n<<<1, 256>>>(N, first, incx, array, incy);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_copy_n(int N, std::complex<double> const* first, int incx, std::complex<double>* array, int incy)
{
  kernel_uninitialized_copy_n<<<1, 256>>>(N, first, incx, array, incy);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_copy_n(int N, int const* first, int incx, int* array, int incy)
{
  kernel_uninitialized_copy_n<<<1, 256>>>(N, first, incx, array, incy);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

// long
void uninitialized_copy_n(long N, double const* first, long incx, double* array, long incy)
{
  kernel_uninitialized_copy_n<<<1, 256>>>(N, first, incx, array, incy);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_copy_n(long N, std::complex<double> const* first, long incx, std::complex<double>* array, long incy)
{
  kernel_uninitialized_copy_n<<<1, 256>>>(N, first, incx, array, incy);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_copy_n(long N, int const* first, long incx, int* array, long incy)
{
  kernel_uninitialized_copy_n<<<1, 256>>>(N, first, incx, array, incy);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

} // namespace kernels
