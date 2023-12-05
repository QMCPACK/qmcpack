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

#define BOOST_NO_AUTO_PTR
#include <cassert>
#include <complex>
#include <type_traits>
#include "AFQMC/Numerics/detail/HIP/hip_kernel_utils.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/hip_settings.h"

namespace kernels
{
template<typename Size, typename Size1, typename T>
__global__ void kernel_fill_n(Size N, T* x, Size1 incx, T const a)
{
  Size N0(8 * blockDim.x * blockIdx.x);
  T* x_(x + Size(incx) * N0);
  Size N_(min(Size(8 * blockDim.x), N - N0));
  for (Size ip = Size(threadIdx.x); ip < N_; ip += Size(blockDim.x))
  {
    x_[ip * Size(incx)] = a;
  }
}

template<typename T, typename Size>
__global__ void kernel_fill2D_n(Size N, Size M, T* y, Size lda, T const a)
{
  for (Size ip = Size(threadIdx.x); ip < N; ip += Size(blockDim.x))
    for (Size jp = Size(threadIdx.y); jp < M; jp += Size(blockDim.y))
    {
      y[ip * lda + jp] = a;
    }
}

void fill_n(char* first, int N, int incx, char const value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  hipLaunchKernelGGL(kernel_fill_n, dim3(nblk), dim3(nthr), 0, 0, N, first, incx, value);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void fill_n(int* first, int N, int incx, int const value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  hipLaunchKernelGGL(kernel_fill_n, dim3(nblk), dim3(nthr), 0, 0, N, first, incx, value);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void fill_n(float* first, int N, int incx, float const value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  hipLaunchKernelGGL(kernel_fill_n, dim3(nblk), dim3(nthr), 0, 0, N, first, incx, value);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void fill_n(double* first, int N, int incx, double const value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  hipLaunchKernelGGL(kernel_fill_n, dim3(nblk), dim3(nthr), 0, 0, N, first, incx, value);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void fill_n(std::complex<float>* first, int N, int incx, std::complex<float> const value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  hipLaunchKernelGGL(kernel_fill_n, dim3(nblk), dim3(nthr), 0, 0, N, first, incx, value);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void fill_n(std::complex<double>* first, int N, int incx, std::complex<double> const value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  hipLaunchKernelGGL(kernel_fill_n, dim3(nblk), dim3(nthr), 0, 0, N, first, incx, value);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void fill_n(char* first, int N, char const value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  hipLaunchKernelGGL(kernel_fill_n, dim3(nblk), dim3(nthr), 0, 0, N, first, 1, value);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void fill_n(long int* first, long unsigned int N, const long int value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  hipLaunchKernelGGL(kernel_fill_n, dim3(nblk), dim3(nthr), 0, 0, N, first, 1, value);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void fill_n(long unsigned int* first, long unsigned int N, const long unsigned int value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  hipLaunchKernelGGL(kernel_fill_n, dim3(nblk), dim3(nthr), 0, 0, N, first, 1, value);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void fill_n(int* first, int N, int const value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  hipLaunchKernelGGL(kernel_fill_n, dim3(nblk), dim3(nthr), 0, 0, N, first, 1, value);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void fill_n(float* first, int N, float const value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  hipLaunchKernelGGL(kernel_fill_n, dim3(nblk), dim3(nthr), 0, 0, N, first, 1, value);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void fill_n(double* first, int N, double const value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  hipLaunchKernelGGL(kernel_fill_n, dim3(nblk), dim3(nthr), 0, 0, N, first, 1, value);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void fill_n(std::complex<float>* first, int N, std::complex<float> const value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  hipLaunchKernelGGL(kernel_fill_n, dim3(nblk), dim3(nthr), 0, 0, N, first, 1, value);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void fill_n(std::complex<double>* first, int N, std::complex<double> const value)
{
  int N_(8 * DEFAULT_BLOCK_SIZE);
  size_t nblk((N + N_ - 1) / N_);
  size_t nthr(DEFAULT_BLOCK_SIZE);
  hipLaunchKernelGGL(kernel_fill_n, dim3(nblk), dim3(nthr), 0, 0, N, first, 1, value);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void fill2D_n(int N, int M, int* A, int lda, int const value)
{
  hipLaunchKernelGGL(kernel_fill2D_n, dim3(32), dim3(32), 0, 0, N, M, A, lda, value);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void fill2D_n(int N, int M, float* A, int lda, float const value)
{
  hipLaunchKernelGGL(kernel_fill2D_n, dim3(32), dim3(32), 0, 0, N, M, A, lda, value);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void fill2D_n(int N, int M, double* A, int lda, double const value)
{
  hipLaunchKernelGGL(kernel_fill2D_n, dim3(32), dim3(32), 0, 0, N, M, A, lda, value);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void fill2D_n(int N, int M, std::complex<double>* A, int lda, std::complex<double> const value)
{
  hipLaunchKernelGGL(kernel_fill2D_n, dim3(32), dim3(32), 0, 0, N, M, A, lda, value);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void fill2D_n(int N, int M, std::complex<float>* A, int lda, std::complex<float> const value)
{
  hipLaunchKernelGGL(kernel_fill2D_n, dim3(32), dim3(32), 0, 0, N, M, A, lda, value);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

} // namespace kernels
