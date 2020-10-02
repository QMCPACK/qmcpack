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
__global__ void kernel_acAxpbB(int m,
                               int n,
                               T const alpha,
                               T const* A,
                               int lda,
                               T const* x,
                               int incx,
                               T const beta,
                               T* B,
                               int ldb)
{
  int i = threadIdx.x + blockDim.x * blockIdx.x;
  int j = threadIdx.y + blockDim.y * blockIdx.y;
  if ((i < m) && (j < n))
    B[j * ldb + i] = beta * B[j * ldb + i] + alpha * A[j * lda + i] * x[i * incx];
}

template<typename T>
__global__ void kernel_acAxpbB(int m,
                               int n,
                               thrust::complex<T> const alpha,
                               thrust::complex<T> const* A,
                               int lda,
                               thrust::complex<T> const* x,
                               int incx,
                               thrust::complex<T> const beta,
                               thrust::complex<T>* B,
                               int ldb)
{
  int i = threadIdx.x + blockDim.x * blockIdx.x;
  int j = threadIdx.y + blockDim.y * blockIdx.y;
  if ((i < m) && (j < n))
    B[j * ldb + i] = beta * B[j * ldb + i] + alpha * conj(A[j * lda + i]) * x[i * incx];
}

void acAxpbB(int m,
             int n,
             double const alpha,
             double const* A,
             int lda,
             double const* x,
             int incx,
             double const beta,
             double* B,
             int ldb)
{
  int xblock_dim = 16;
  int xgrid_dim  = (m + xblock_dim - 1) / xblock_dim;
  int yblock_dim = 16;
  int ygrid_dim  = (n + yblock_dim - 1) / yblock_dim;
  dim3 block_dim(xblock_dim, yblock_dim);
  dim3 grid_dim(xgrid_dim, ygrid_dim);
  hipLaunchKernelGGL(kernel_acAxpbB, dim3(grid_dim), dim3(block_dim), 0, 0, m, n, alpha, A, lda, x, incx, beta, B, ldb);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void acAxpbB(int m,
             int n,
             float const alpha,
             float const* A,
             int lda,
             float const* x,
             int incx,
             float const beta,
             float* B,
             int ldb)
{
  int xblock_dim = 16;
  int xgrid_dim  = (m + xblock_dim - 1) / xblock_dim;
  int yblock_dim = 16;
  int ygrid_dim  = (n + yblock_dim - 1) / yblock_dim;
  dim3 block_dim(xblock_dim, yblock_dim);
  dim3 grid_dim(xgrid_dim, ygrid_dim);
  hipLaunchKernelGGL(kernel_acAxpbB, dim3(grid_dim), dim3(block_dim), 0, 0, m, n, alpha, A, lda, x, incx, beta, B, ldb);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void acAxpbB(int m,
             int n,
             std::complex<double> const alpha,
             std::complex<double> const* A,
             int lda,
             std::complex<double> const* x,
             int incx,
             std::complex<double> const beta,
             std::complex<double>* B,
             int ldb)
{
  int xblock_dim = 16;
  int xgrid_dim  = (m + xblock_dim - 1) / xblock_dim;
  int yblock_dim = 16;
  int ygrid_dim  = (n + yblock_dim - 1) / yblock_dim;
  dim3 block_dim(xblock_dim, yblock_dim);
  dim3 grid_dim(xgrid_dim, ygrid_dim);
  hipLaunchKernelGGL(kernel_acAxpbB, dim3(grid_dim), dim3(block_dim), 0, 0, m, n,
                     static_cast<thrust::complex<double> const>(alpha),
                     reinterpret_cast<thrust::complex<double> const*>(A), lda,
                     reinterpret_cast<thrust::complex<double> const*>(x), incx,
                     static_cast<thrust::complex<double> const>(beta), reinterpret_cast<thrust::complex<double>*>(B),
                     ldb);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void acAxpbB(int m,
             int n,
             std::complex<float> const alpha,
             std::complex<float> const* A,
             int lda,
             std::complex<float> const* x,
             int incx,
             std::complex<float> const beta,
             std::complex<float>* B,
             int ldb)
{
  int xblock_dim = 16;
  int xgrid_dim  = (m + xblock_dim - 1) / xblock_dim;
  int yblock_dim = 16;
  int ygrid_dim  = (n + yblock_dim - 1) / yblock_dim;
  dim3 block_dim(xblock_dim, yblock_dim);
  dim3 grid_dim(xgrid_dim, ygrid_dim);
  hipLaunchKernelGGL(kernel_acAxpbB, dim3(grid_dim), dim3(block_dim), 0, 0, m, n,
                     static_cast<thrust::complex<float> const>(alpha),
                     reinterpret_cast<thrust::complex<float> const*>(A), lda,
                     reinterpret_cast<thrust::complex<float> const*>(x), incx,
                     static_cast<thrust::complex<float> const>(beta), reinterpret_cast<thrust::complex<float>*>(B),
                     ldb);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}


} // namespace kernels
