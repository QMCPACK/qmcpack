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
template<typename T, typename Q>
__global__ void kernel_copy_n_cast(T const* A, int n, Q* B)
{
  int i = threadIdx.x + blockDim.x * blockIdx.x;
  if (i < n)
    B[i] = static_cast<Q>(A[i]);
}

template<typename T, typename Q>
__global__ void kernel_copy_n_cast(thrust::complex<T> const* A, int n, thrust::complex<Q>* B)
{
  int i = threadIdx.x + blockDim.x * blockIdx.x;
  if (i < n)
    B[i] = static_cast<thrust::complex<Q>>(A[i]);
}

void copy_n_cast(double const* A, int n, float* B)
{
  int block_dim = 256;
  int grid_dim  = (n + block_dim - 1) / block_dim;
  hipLaunchKernelGGL(kernel_copy_n_cast, dim3(grid_dim), dim3(block_dim), 0, 0, A, n, B);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void copy_n_cast(float const* A, int n, double* B)
{
  int block_dim = 256;
  int grid_dim  = (n + block_dim - 1) / block_dim;
  hipLaunchKernelGGL(kernel_copy_n_cast, dim3(grid_dim), dim3(block_dim), 0, 0, A, n, B);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void copy_n_cast(std::complex<double> const* A, int n, std::complex<float>* B)
{
  int block_dim = 256;
  int grid_dim  = (n + block_dim - 1) / block_dim;
  hipLaunchKernelGGL(kernel_copy_n_cast, dim3(grid_dim), dim3(block_dim), 0, 0,
                     reinterpret_cast<thrust::complex<double> const*>(A), n,
                     reinterpret_cast<thrust::complex<float>*>(B));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void copy_n_cast(std::complex<float> const* A, int n, std::complex<double>* B)
{
  int block_dim = 256;
  int grid_dim  = (n + block_dim - 1) / block_dim;
  hipLaunchKernelGGL(kernel_copy_n_cast, dim3(grid_dim), dim3(block_dim), 0, 0,
                     reinterpret_cast<thrust::complex<float> const*>(A), n,
                     reinterpret_cast<thrust::complex<double>*>(B));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

} // namespace kernels
