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
__global__ void kernel_setIdentity(int m, int n, T* A, int lda)
{
  int i = threadIdx.x + blockDim.x * blockIdx.x;
  int j = threadIdx.y + blockDim.y * blockIdx.y;
  if ((i < m) && (j < n))
  {
    if (i == j)
    {
      A[i * lda + i] = T(1.0);
    }
    else
    {
      A[j * lda + i] = T(0.0);
    }
  }
}

template<typename T>
__global__ void kernel_setIdentity(int m, int n, thrust::complex<T>* A, int lda)
{
  int i = threadIdx.x + blockDim.x * blockIdx.x;
  int j = threadIdx.y + blockDim.y * blockIdx.y;
  if ((i < m) && (j < n))
  {
    if (i == j)
    {
      A[i * lda + i] = thrust::complex<T>(1.0, 0.0);
    }
    else
    {
      A[j * lda + i] = thrust::complex<T>(0.0, 0.0);
    }
  }
}

template<typename T>
__global__ void kernel_setIdentity_strided(int nbatch, int stride, int m, int n, T* A, int lda)
{
  int i = threadIdx.x + blockDim.x * blockIdx.x;
  int j = threadIdx.y + blockDim.y * blockIdx.y;
  if ((i < m) && (j < n) && (blockIdx.z < nbatch))
  {
    if (i == j)
    {
      A[blockIdx.z * stride + i * lda + i] = T(1.0);
    }
    else
    {
      A[blockIdx.z * stride + j * lda + i] = T(0.0);
    }
  }
}

template<typename T>
__global__ void kernel_setIdentity_strided(int nbatch, int stride, int m, int n, thrust::complex<T>* A, int lda)
{
  int i = threadIdx.x + blockDim.x * blockIdx.x;
  int j = threadIdx.y + blockDim.y * blockIdx.y;
  if ((i < m) && (j < n) && (blockIdx.z < nbatch))
  {
    if (i == j)
    {
      A[blockIdx.z * stride + i * lda + i] = thrust::complex<T>(1.0, 0.0);
    }
    else
    {
      A[blockIdx.z * stride + j * lda + i] = thrust::complex<T>(0.0, 0.0);
    }
  }
}

void set_identity(int m, int n, double* A, int lda)
{
  int xblock_dim = 16;
  int xgrid_dim  = (m + xblock_dim - 1) / xblock_dim;
  int ygrid_dim  = (n + xblock_dim - 1) / xblock_dim;
  dim3 block_dim(xblock_dim, xblock_dim);
  dim3 grid_dim(xgrid_dim, ygrid_dim);
  hipLaunchKernelGGL(kernel_setIdentity, dim3(grid_dim), dim3(block_dim), 0, 0, m, n, A, lda);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void set_identity(int m, int n, float* A, int lda)
{
  int xblock_dim = 16;
  int xgrid_dim  = (m + xblock_dim - 1) / xblock_dim;
  int ygrid_dim  = (n + xblock_dim - 1) / xblock_dim;
  dim3 block_dim(xblock_dim, xblock_dim);
  dim3 grid_dim(xgrid_dim, ygrid_dim);
  hipLaunchKernelGGL(kernel_setIdentity, dim3(grid_dim), dim3(block_dim), 0, 0, m, n, A, lda);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void set_identity(int m, int n, std::complex<double>* A, int lda)
{
  int xblock_dim = 16;
  int xgrid_dim  = (m + xblock_dim - 1) / xblock_dim;
  int ygrid_dim  = (n + xblock_dim - 1) / xblock_dim;
  dim3 block_dim(xblock_dim, xblock_dim);
  dim3 grid_dim(xgrid_dim, ygrid_dim);
  hipLaunchKernelGGL(kernel_setIdentity, dim3(grid_dim), dim3(block_dim), 0, 0, m, n,
                     reinterpret_cast<thrust::complex<double>*>(A), lda);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void set_identity(int m, int n, std::complex<float>* A, int lda)
{
  int xblock_dim = 16;
  int xgrid_dim  = (m + xblock_dim - 1) / xblock_dim;
  int ygrid_dim  = (n + xblock_dim - 1) / xblock_dim;
  dim3 block_dim(xblock_dim, xblock_dim);
  dim3 grid_dim(xgrid_dim, ygrid_dim);
  hipLaunchKernelGGL(kernel_setIdentity, dim3(grid_dim), dim3(block_dim), 0, 0, m, n,
                     reinterpret_cast<thrust::complex<float>*>(A), lda);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void set_identity_strided(int nbatch, int stride, int m, int n, double* A, int lda)
{
  int xblock_dim = 16;
  int xgrid_dim  = (m + xblock_dim - 1) / xblock_dim;
  int ygrid_dim  = (n + xblock_dim - 1) / xblock_dim;
  dim3 block_dim(xblock_dim, xblock_dim);
  dim3 grid_dim(xgrid_dim, ygrid_dim, nbatch);
  hipLaunchKernelGGL(kernel_setIdentity_strided, dim3(grid_dim), dim3(block_dim), 0, 0, nbatch, stride, m, n, A, lda);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void set_identity_strided(int nbatch, int stride, int m, int n, float* A, int lda)
{
  int xblock_dim = 16;
  int xgrid_dim  = (m + xblock_dim - 1) / xblock_dim;
  int ygrid_dim  = (n + xblock_dim - 1) / xblock_dim;
  dim3 block_dim(xblock_dim, xblock_dim);
  dim3 grid_dim(xgrid_dim, ygrid_dim, nbatch);
  hipLaunchKernelGGL(kernel_setIdentity_strided, dim3(grid_dim), dim3(block_dim), 0, 0, nbatch, stride, m, n, A, lda);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void set_identity_strided(int nbatch, int stride, int m, int n, std::complex<double>* A, int lda)
{
  int xblock_dim = 16;
  int xgrid_dim  = (m + xblock_dim - 1) / xblock_dim;
  int ygrid_dim  = (n + xblock_dim - 1) / xblock_dim;
  dim3 block_dim(xblock_dim, xblock_dim);
  dim3 grid_dim(xgrid_dim, ygrid_dim, nbatch);
  hipLaunchKernelGGL(kernel_setIdentity_strided, dim3(grid_dim), dim3(block_dim), 0, 0, nbatch, stride, m, n,
                     reinterpret_cast<thrust::complex<double>*>(A), lda);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void set_identity_strided(int nbatch, int stride, int m, int n, std::complex<float>* A, int lda)
{
  int xblock_dim = 16;
  int xgrid_dim  = (m + xblock_dim - 1) / xblock_dim;
  int ygrid_dim  = (n + xblock_dim - 1) / xblock_dim;
  dim3 block_dim(xblock_dim, xblock_dim);
  dim3 grid_dim(xgrid_dim, ygrid_dim, nbatch);
  hipLaunchKernelGGL(kernel_setIdentity_strided, dim3(grid_dim), dim3(block_dim), 0, 0, nbatch, stride, m, n,
                     reinterpret_cast<thrust::complex<float>*>(A), lda);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

} // namespace kernels
