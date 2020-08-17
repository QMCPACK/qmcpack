///////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Fionn Malone, malone14@llnl.gov, LLNL
//
// File created by: Fionn Malone, malone14@llnl.gov, LLNL
////////////////////////////////////////////////////////////////////////////////


#include <cassert>
#include <complex>
#include <type_traits>
#include <thrust/complex.h>
#include "AFQMC/Memory/HIP/hip_utilities.h"
//#include "AFQMC/Numerics/detail/HIP/Kernels/strided_range.hpp"

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
  hipLaunchKernelGGL(kernel_uninitialized_copy_n, dim3(1), dim3(256), 0, 0, N, first, incx, array, incy);
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

void uninitialized_copy_n(int N, std::complex<double> const* first, int incx, std::complex<double>* array, int incy)
{
  hipLaunchKernelGGL(kernel_uninitialized_copy_n, dim3(1), dim3(256), 0, 0, N, first, incx, array, incy);
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

void uninitialized_copy_n(int N, int const* first, int incx, int* array, int incy)
{
  hipLaunchKernelGGL(kernel_uninitialized_copy_n, dim3(1), dim3(256), 0, 0, N, first, incx, array, incy);
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

// long
void uninitialized_copy_n(long N, double const* first, long incx, double* array, long incy)
{
  hipLaunchKernelGGL(kernel_uninitialized_copy_n, dim3(1), dim3(256), 0, 0, N, first, incx, array, incy);
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

void uninitialized_copy_n(long N, std::complex<double> const* first, long incx, std::complex<double>* array, long incy)
{
  hipLaunchKernelGGL(kernel_uninitialized_copy_n, dim3(1), dim3(256), 0, 0, N, first, incx, array, incy);
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

void uninitialized_copy_n(long N, int const* first, long incx, int* array, long incy)
{
  hipLaunchKernelGGL(kernel_uninitialized_copy_n, dim3(1), dim3(256), 0, 0, N, first, incx, array, incy);
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

} // namespace kernels
