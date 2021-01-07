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
// Tab [nbatch][nwalk][nocc][nocc][nchol]
template<typename T, typename T2>
__global__ void kernel_batched_dot_wabn_wban(int nbatch,
                                             int nwalk,
                                             int nocc,
                                             int nchol,
                                             thrust::complex<T2> const* alpha,
                                             thrust::complex<T2> const* Tab,
                                             thrust::complex<T>* y,
                                             int incy)
{
  int batch = blockIdx.x;
  if (batch >= nbatch)
    return;
  if (blockIdx.y >= nwalk * nocc * nocc)
    return;
  __shared__ uninitialized_array<thrust::complex<T>, DOT_BLOCK_SIZE> cache;
  int nocc2              = nocc * nocc;
  int w                  = blockIdx.y / (nocc2);
  int a                  = (blockIdx.y % (nocc2)) / nocc;
  int b                  = (blockIdx.y % (nocc2)) % nocc;
  int i                  = threadIdx.x;
  thrust::complex<T> alp = static_cast<thrust::complex<T>>(alpha[batch]);
  thrust::complex<T2> const* A_(Tab + 2 * batch * nwalk * nocc2 * nchol + ((w * nocc + a) * nocc + b) * nchol);
  thrust::complex<T2> const* B_(Tab + (2 * batch + 1) * nwalk * nocc2 * nchol + ((w * nocc + b) * nocc + a) * nchol);
  cache[threadIdx.x] = thrust::complex<T>(0.0);
  while (i < nchol)
  {
    cache[threadIdx.x] += static_cast<thrust::complex<T>>(A_[i] * B_[i]);
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
  //if( threadIdx.x == 0 ) *(y+w*incy) = (*(y+w*incy)) + alp * cache[ 0 ];
  if (threadIdx.x == 0)
  {
    T re   = (alp * cache[0]).real();
    T im   = (alp * cache[0]).imag();
    T* re_ = reinterpret_cast<T*>(y + w * incy);
    atomicAdd(re_, re);
    atomicAdd(re_ + 1, im);
  }
}

template<typename T, typename T2>
__global__ void kernel_batched_dot_wanb_wbna(int nbatch,
                                             int nwalk,
                                             int nocc,
                                             int nchol,
                                             thrust::complex<T2> const* alpha,
                                             thrust::complex<T2> const* Tab,
                                             thrust::complex<T>* y,
                                             int incy)
{
  int batch = blockIdx.x;
  if (batch >= nbatch)
    return;
  if (blockIdx.y >= nwalk * nocc * nocc)
    return;
  __shared__ uninitialized_array<thrust::complex<T>, DOT_BLOCK_SIZE> cache;
  int nocc2              = nocc * nocc;
  int w                  = blockIdx.y / (nocc2);
  int a                  = (blockIdx.y % (nocc2)) / nocc;
  int b                  = (blockIdx.y % (nocc2)) % nocc;
  int i                  = threadIdx.x;
  thrust::complex<T> alp = static_cast<thrust::complex<T>>(alpha[batch]);
  thrust::complex<T2> const* A_(Tab + 2 * batch * nwalk * nocc2 * nchol + ((w * nocc + a) * nocc) * nchol + b);
  thrust::complex<T2> const* B_(Tab + (2 * batch + 1) * nwalk * nocc2 * nchol + ((w * nocc + b) * nocc) * nchol + a);
  cache[threadIdx.x] = thrust::complex<T>(0.0);
  while (i < nchol)
  {
    cache[threadIdx.x] += static_cast<thrust::complex<T>>(A_[i * nocc] * B_[i * nocc]);
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
  //if( threadIdx.x == 0 ) *(y+w*incy) = (*(y+w*incy)) + alp * cache[ 0 ];
  if (threadIdx.x == 0)
  {
    T re   = (alp * cache[0]).real();
    T im   = (alp * cache[0]).imag();
    T* re_ = reinterpret_cast<T*>(y + w * incy);
    atomicAdd(re_, re);
    atomicAdd(re_ + 1, im);
  }
}

void batched_dot_wabn_wban(int nbatch,
                           int nwalk,
                           int nocc,
                           int nchol,
                           std::complex<double> const* alpha,
                           std::complex<double> const* Tab,
                           std::complex<double>* y,
                           int incy)
{
  int n_ = nwalk * nocc * nocc;
  dim3 grid_dim(nbatch, n_, 1);
  hipLaunchKernelGGL(kernel_batched_dot_wabn_wban, dim3(grid_dim), dim3(DOT_BLOCK_SIZE), 0, 0, nbatch, nwalk, nocc,
                     nchol, reinterpret_cast<thrust::complex<double> const*>(alpha),
                     reinterpret_cast<thrust::complex<double> const*>(Tab),
                     reinterpret_cast<thrust::complex<double>*>(y), incy);
  qmc_hip::hip_kernel_check(hipGetLastError(), "batched_dot_wabn_wban");
  qmc_hip::hip_kernel_check(hipDeviceSynchronize(), "batched_dot_wabn_wban");
}

void batched_dot_wabn_wban(int nbatch,
                           int nwalk,
                           int nocc,
                           int nchol,
                           std::complex<float> const* alpha,
                           std::complex<float> const* Tab,
                           std::complex<float>* y,
                           int incy)
{
  int n_ = nwalk * nocc * nocc;
  dim3 grid_dim(nbatch, n_, 1);
  hipLaunchKernelGGL(kernel_batched_dot_wabn_wban, dim3(grid_dim), dim3(DOT_BLOCK_SIZE), 0, 0, nbatch, nwalk, nocc,
                     nchol, reinterpret_cast<thrust::complex<float> const*>(alpha),
                     reinterpret_cast<thrust::complex<float> const*>(Tab), reinterpret_cast<thrust::complex<float>*>(y),
                     incy);
  qmc_hip::hip_kernel_check(hipGetLastError(), "batched_dot_wabn_wban");
  qmc_hip::hip_kernel_check(hipDeviceSynchronize(), "batched_dot_wabn_wban");
}

void batched_dot_wabn_wban(int nbatch,
                           int nwalk,
                           int nocc,
                           int nchol,
                           std::complex<float> const* alpha,
                           std::complex<float> const* Tab,
                           std::complex<double>* y,
                           int incy)
{
  int n_ = nwalk * nocc * nocc;
  dim3 grid_dim(nbatch, n_, 1);
  hipLaunchKernelGGL(kernel_batched_dot_wabn_wban, dim3(grid_dim), dim3(DOT_BLOCK_SIZE), 0, 0, nbatch, nwalk, nocc,
                     nchol, reinterpret_cast<thrust::complex<float> const*>(alpha),
                     reinterpret_cast<thrust::complex<float> const*>(Tab),
                     reinterpret_cast<thrust::complex<double>*>(y), incy);
  qmc_hip::hip_kernel_check(hipGetLastError(), "batched_dot_wabn_wban");
  qmc_hip::hip_kernel_check(hipDeviceSynchronize(), "batched_dot_wabn_wban");
}

// anb/bna
void batched_dot_wanb_wbna(int nbatch,
                           int nwalk,
                           int nocc,
                           int nchol,
                           std::complex<double> const* alpha,
                           std::complex<double> const* Tab,
                           std::complex<double>* y,
                           int incy)
{
  int n_ = nwalk * nocc * nocc;
  dim3 grid_dim(nbatch, n_, 1);
  hipLaunchKernelGGL(kernel_batched_dot_wanb_wbna, dim3(grid_dim), dim3(DOT_BLOCK_SIZE), 0, 0, nbatch, nwalk, nocc,
                     nchol, reinterpret_cast<thrust::complex<double> const*>(alpha),
                     reinterpret_cast<thrust::complex<double> const*>(Tab),
                     reinterpret_cast<thrust::complex<double>*>(y), incy);
  qmc_hip::hip_kernel_check(hipGetLastError(), "batched_dot_wanb_wbna");
  qmc_hip::hip_kernel_check(hipDeviceSynchronize(), "batched_dot_wanb_wbna");
}

void batched_dot_wanb_wbna(int nbatch,
                           int nwalk,
                           int nocc,
                           int nchol,
                           std::complex<float> const* alpha,
                           std::complex<float> const* Tab,
                           std::complex<float>* y,
                           int incy)
{
  int n_ = nwalk * nocc * nocc;
  dim3 grid_dim(nbatch, n_, 1);
  hipLaunchKernelGGL(kernel_batched_dot_wanb_wbna, dim3(grid_dim), dim3(DOT_BLOCK_SIZE), 0, 0, nbatch, nwalk, nocc,
                     nchol, reinterpret_cast<thrust::complex<float> const*>(alpha),
                     reinterpret_cast<thrust::complex<float> const*>(Tab), reinterpret_cast<thrust::complex<float>*>(y),
                     incy);
  qmc_hip::hip_kernel_check(hipGetLastError(), "batched_dot_wanb_wbna");
  qmc_hip::hip_kernel_check(hipDeviceSynchronize(), "batched_dot_wanb_wbna");
}

void batched_dot_wanb_wbna(int nbatch,
                           int nwalk,
                           int nocc,
                           int nchol,
                           std::complex<float> const* alpha,
                           std::complex<float> const* Tab,
                           std::complex<double>* y,
                           int incy)
{
  int n_ = nwalk * nocc * nocc;
  dim3 grid_dim(nbatch, n_, 1);
  hipLaunchKernelGGL(kernel_batched_dot_wanb_wbna, dim3(grid_dim), dim3(DOT_BLOCK_SIZE), 0, 0, nbatch, nwalk, nocc,
                     nchol, reinterpret_cast<thrust::complex<float> const*>(alpha),
                     reinterpret_cast<thrust::complex<float> const*>(Tab),
                     reinterpret_cast<thrust::complex<double>*>(y), incy);
  qmc_hip::hip_kernel_check(hipGetLastError(), "batched_dot_wanb_wbna");
  qmc_hip::hip_kernel_check(hipDeviceSynchronize(), "batched_dot_wanb_wbna");
}


} // namespace kernels
