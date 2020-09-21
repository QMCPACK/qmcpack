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
#include "uninitialized_array.hpp"
#include "AFQMC/Numerics/detail/HIP/Kernels/hip_settings.h"
#include "AFQMC/Numerics/detail/HIP/hip_kernel_utils.h"

namespace kernels
{
template<typename T>
__global__ void kernel_dot(int n, T const alpha, T const* A, int lda, T const* B, int ldb, T const beta, T* y, int incy)
{
  __shared__ uninitialized_array<T, DOT_BLOCK_SIZE> cache;
  int i              = threadIdx.x;
  int j              = blockIdx.x;
  cache[threadIdx.x] = T(0.0);
  while (i < n)
  {
    cache[threadIdx.x] += A[j * lda + i] * B[j * ldb + i];
    i += blockDim.x;
  }
  __syncthreads(); // required because later on the hiprrent thread is accessing
                   // data written by another thread
  i = DOT_BLOCK_SIZE / 2;
  while (i > 0)
  {
    if (threadIdx.x < i)
      cache[threadIdx.x] += cache[threadIdx.x + i];
    __syncthreads();
    i /= 2; //not sure bitwise operations are actually faster
  }
  if (threadIdx.x == 0)
    *(y + j * incy) = beta * (*(y + j * incy)) + alpha * cache[0];
}

template<typename T>
__global__ void kernel_dot(int n,
                           thrust::complex<T> const alpha,
                           thrust::complex<T> const* A,
                           int lda,
                           thrust::complex<T> const* B,
                           int ldb,
                           thrust::complex<T> const beta,
                           thrust::complex<T>* y,
                           int incy)
{
  __shared__ uninitialized_array<thrust::complex<T>, DOT_BLOCK_SIZE> cache;
  int i              = threadIdx.x;
  int j              = blockIdx.x;
  cache[threadIdx.x] = thrust::complex<T>(0.0, 0.0);
  while (i < n)
  {
    cache[threadIdx.x] += A[j * lda + i] * B[j * ldb + i];
    i += blockDim.x;
  }
  __syncthreads(); // required because later on the hiprrent thread is accessing
                   // data written by another thread
  i = DOT_BLOCK_SIZE / 2;
  while (i > 0)
  {
    if (threadIdx.x < i)
      cache[threadIdx.x] += cache[threadIdx.x + i];
    __syncthreads();
    i /= 2; //not sure bitwise operations are actually faster
  }
  if (threadIdx.x == 0)
    *(y + j * incy) = beta * (*(y + j * incy)) + alpha * cache[0];
}

// y[i] = beta * y[i] + sum_k alpha * A[k,i] * B[k,i]
void batchedDot(int m,
                int n,
                double const alpha,
                double const* A,
                int lda,
                double const* B,
                int ldb,
                double const beta,
                double* y,
                int incy)
{
  hipLaunchKernelGGL(kernel_dot, dim3(m), dim3(DOT_BLOCK_SIZE), 0, 0, n, alpha, A, lda, B, ldb, beta, y, incy);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void batchedDot(int m,
                int n,
                float const alpha,
                float const* A,
                int lda,
                float const* B,
                int ldb,
                float const beta,
                float* y,
                int incy)
{
  hipLaunchKernelGGL(kernel_dot, dim3(m), dim3(DOT_BLOCK_SIZE), 0, 0, n, alpha, A, lda, B, ldb, beta, y, incy);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void batchedDot(int m,
                int n,
                std::complex<double> const alpha,
                std::complex<double> const* A,
                int lda,
                std::complex<double> const* B,
                int ldb,
                std::complex<double> const beta,
                std::complex<double>* y,
                int incy)
{
  hipLaunchKernelGGL(kernel_dot, dim3(m), dim3(DOT_BLOCK_SIZE), 0, 0, n,
                     static_cast<thrust::complex<double> const>(alpha),
                     reinterpret_cast<thrust::complex<double> const*>(A), lda,
                     reinterpret_cast<thrust::complex<double> const*>(B), ldb,
                     static_cast<thrust::complex<double> const>(beta), reinterpret_cast<thrust::complex<double>*>(y),
                     incy);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void batchedDot(int m,
                int n,
                std::complex<float> const alpha,
                std::complex<float> const* A,
                int lda,
                std::complex<float> const* B,
                int ldb,
                std::complex<float> const beta,
                std::complex<float>* y,
                int incy)
{
  hipLaunchKernelGGL(kernel_dot, dim3(m), dim3(DOT_BLOCK_SIZE), 0, 0, n,
                     static_cast<thrust::complex<float> const>(alpha),
                     reinterpret_cast<thrust::complex<float> const*>(A), lda,
                     reinterpret_cast<thrust::complex<float> const*>(B), ldb,
                     static_cast<thrust::complex<float> const>(beta), reinterpret_cast<thrust::complex<float>*>(y),
                     incy);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

} // namespace kernels
