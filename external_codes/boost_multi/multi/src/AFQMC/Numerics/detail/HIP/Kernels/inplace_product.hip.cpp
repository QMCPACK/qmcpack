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
#include <thrust/complex.h>
#include <hip/hip_runtime.h>
#include "AFQMC/Numerics/detail/HIP/hip_kernel_utils.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/hip_settings.h"

namespace kernels
{
template<typename T1, typename T>
__global__ void kernel_inplace_product(int nbatch, int n, int m, T1 const* B, int ldb, T* A, int lda)
{
  int nb = blockIdx.z;
  int i  = blockIdx.x;
  if ((i < n) && (nb < nbatch))
  {
    int j = blockDim.y * blockIdx.y + threadIdx.y;
    if (j < m)
      *(A + (nb * n + i) * lda + j) *= static_cast<T>((*(B + i * ldb + j)));
  }
}

void inplace_product(int nbatch, int n, int m, double const* B, int ldb, std::complex<double>* A, int lda)
{
  int nby = (m + MAX_THREADS_PER_DIM - 1L) / MAX_THREADS_PER_DIM;
  dim3 grid_dim(n, nby, nbatch);
  dim3 block_dim(1, MAX_THREADS_PER_DIM, 1);
  hipLaunchKernelGGL(kernel_inplace_product, dim3(grid_dim), dim3(block_dim), 0, 0, nbatch, n, m, B, ldb,
                     reinterpret_cast<thrust::complex<double>*>(A), lda);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void inplace_product(int nbatch, int n, int m, std::complex<double> const* B, int ldb, std::complex<double>* A, int lda)
{
  int nby = (m + MAX_THREADS_PER_DIM - 1L) / MAX_THREADS_PER_DIM;
  dim3 grid_dim(n, nby, nbatch);
  dim3 block_dim(1, MAX_THREADS_PER_DIM, 1);
  hipLaunchKernelGGL(kernel_inplace_product, dim3(grid_dim), dim3(block_dim), 0, 0, nbatch, n, m,
                     reinterpret_cast<thrust::complex<double> const*>(B), ldb,
                     reinterpret_cast<thrust::complex<double>*>(A), lda);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void inplace_product(int nbatch, int n, int m, float const* B, int ldb, std::complex<float>* A, int lda)
{
  int nby = (m + MAX_THREADS_PER_DIM - 1L) / MAX_THREADS_PER_DIM;
  dim3 grid_dim(n, nby, nbatch);
  dim3 block_dim(1, MAX_THREADS_PER_DIM, 1);
  hipLaunchKernelGGL(kernel_inplace_product, dim3(grid_dim), dim3(block_dim), 0, 0, nbatch, n, m, B, ldb,
                     reinterpret_cast<thrust::complex<float>*>(A), lda);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void inplace_product(int nbatch, int n, int m, std::complex<float> const* B, int ldb, std::complex<float>* A, int lda)
{
  int nby = (m + MAX_THREADS_PER_DIM - 1L) / MAX_THREADS_PER_DIM;
  dim3 grid_dim(n, nby, nbatch);
  dim3 block_dim(1, MAX_THREADS_PER_DIM, 1);
  hipLaunchKernelGGL(kernel_inplace_product, dim3(grid_dim), dim3(block_dim), 0, 0, nbatch, n, m,
                     reinterpret_cast<thrust::complex<float> const*>(B), ldb,
                     reinterpret_cast<thrust::complex<float>*>(A), lda);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

} // namespace kernels
