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
#include <hip/hip_runtime.h>
#include <thrust/complex.h>
#include <hip/hip_runtime.h>
#include "AFQMC/Numerics/detail/HIP/hip_kernel_utils.h"

namespace kernels
{
template<typename T>
__global__ void kernel_adiagApy(int N, T const alpha, T const* A, int lda, T* y, int incy)
{
  int i = threadIdx.x + blockDim.x * blockIdx.x;
  if (i < N)
  {
    y[i * incy] += alpha * A[i * lda + i];
  }
}

template<typename T>
__global__ void kernel_adiagApy(int N,
                                thrust::complex<T> const alpha,
                                thrust::complex<T> const* A,
                                int lda,
                                thrust::complex<T>* y,
                                int incy)
{
  int i = threadIdx.x + blockDim.x * blockIdx.x;
  if (i < N)
  {
    y[i * incy] += alpha * A[i * lda + i];
  }
}

void adiagApy(int N, double const alpha, double const* A, int lda, double* y, int incy)
{
  int block_dim = 256;
  int grid_dim  = (N + block_dim - 1) / block_dim;
  hipLaunchKernelGGL(kernel_adiagApy, dim3(grid_dim), dim3(block_dim), 0, 0, N, alpha, A, lda, y, incy);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void adiagApy(int N,
              std::complex<double> const alpha,
              std::complex<double> const* A,
              int lda,
              std::complex<double>* y,
              int incy)
{
  int block_dim = 256;
  int grid_dim  = (N + block_dim - 1) / block_dim;
  hipLaunchKernelGGL(kernel_adiagApy, dim3(grid_dim), dim3(block_dim), 0, 0, N,
                     static_cast<thrust::complex<double> const>(alpha),
                     reinterpret_cast<thrust::complex<double> const*>(A), lda,
                     reinterpret_cast<thrust::complex<double>*>(y), incy);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void adiagApy(int N, float const alpha, float const* A, int lda, float* y, int incy)
{
  int block_dim = 256;
  int grid_dim  = (N + block_dim - 1) / block_dim;
  hipLaunchKernelGGL(kernel_adiagApy, dim3(grid_dim), dim3(block_dim), 0, 0, N, alpha, A, lda, y, incy);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void adiagApy(int N,
              std::complex<float> const alpha,
              std::complex<float> const* A,
              int lda,
              std::complex<float>* y,
              int incy)
{
  int block_dim = 128;
  int grid_dim  = (N + block_dim - 1) / block_dim;
  hipLaunchKernelGGL(kernel_adiagApy, dim3(grid_dim), dim3(block_dim), 0, 0, N,
                     static_cast<thrust::complex<float> const>(alpha),
                     reinterpret_cast<thrust::complex<float> const*>(A), lda,
                     reinterpret_cast<thrust::complex<float>*>(y), incy);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

} // namespace kernels
