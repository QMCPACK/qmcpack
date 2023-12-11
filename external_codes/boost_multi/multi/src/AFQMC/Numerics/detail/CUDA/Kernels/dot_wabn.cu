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
#include <cuda.h>
#include <thrust/complex.h>
#include <cuda_runtime.h>
#include <thrust/system/cuda/detail/core/util.h>
#include "AFQMC/Numerics/detail/CUDA/Kernels/cuda_settings.h"
#define ENABLE_CUDA 1
#include "AFQMC/Memory/CUDA/cuda_utilities.h"
#if __CUDA_ARCH__ < 600
#include "AFQMC/Numerics/detail/CUDA/Kernels/cuda_workaround_legacy_hardware.cuh"
#endif

namespace kernels
{
// Tab nwalk][nocc][nocc][nchol]
template<typename T, typename T2>
__global__ void kernel_dot_wabn(int nwalk,
                                int nocc,
                                int nchol,
                                thrust::complex<T2> const alpha,
                                thrust::complex<T2> const* Tab,
                                thrust::complex<T>* y,
                                int incy)
{
  if (blockIdx.x >= nwalk * nocc * nocc)
    return;
  __shared__ thrust::cuda_cub::core::uninitialized_array<thrust::complex<T>, DOT_BLOCK_SIZE> cache;
  int nocc2              = nocc * nocc;
  int w                  = blockIdx.x / (nocc2);
  int a                  = (blockIdx.x % (nocc2)) / nocc;
  int b                  = (blockIdx.x % (nocc2)) % nocc;
  int i                  = threadIdx.x;
  thrust::complex<T> alp = static_cast<thrust::complex<T>>(alpha);
  thrust::complex<T2> const* A_(Tab + ((w * nocc + a) * nocc + b) * nchol);
  thrust::complex<T2> const* B_(Tab + ((w * nocc + b) * nocc + a) * nchol);
  cache[threadIdx.x] = thrust::complex<T>(0.0);
  while (i < nchol)
  {
    cache[threadIdx.x] += static_cast<thrust::complex<T>>(A_[i] * B_[i]);
    i += blockDim.x;
  }
  __syncthreads(); // required because later on the current thread is accessing
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
  {
    T re   = (alp * cache[0]).real();
    T im   = (alp * cache[0]).imag();
    T* re_ = reinterpret_cast<T*>(y + w * incy);
    atomicAdd(re_, re);
    atomicAdd(re_ + 1, im);
  }
}

template<typename T, typename T2>
__global__ void kernel_dot_wanb(int nt,
                                int nwalk,
                                int nocc,
                                int nchol,
                                thrust::complex<T2> const alpha,
                                thrust::complex<T2> const* Tab,
                                thrust::complex<T>* y,
                                int incy)
{
  if (blockIdx.x >= nwalk)
    return;

  int a  = blockIdx.y * blockDim.x + threadIdx.x;
  int nb = blockIdx.z * blockDim.y * nt + threadIdx.y;

  __shared__ thrust::cuda_cub::core::uninitialized_array<thrust::complex<T>, 1024> cache;

  int nid = blockDim.x * blockDim.y * blockDim.z;
  int id  = (threadIdx.x * blockDim.y + threadIdx.y) * blockDim.z + threadIdx.z;

  cache[id]              = thrust::complex<T>(0.0);
  thrust::complex<T> alp = static_cast<thrust::complex<T>>(alpha);
  thrust::complex<T2> const* A_(Tab + blockIdx.x * nocc * nocc * nchol);
  for (int n = 0; n < nt; n++, nb += blockDim.y)
  {
    if (a >= nocc || nb >= nocc * nchol)
      break;
    int i = nb / nocc;
    int b = nb % nocc;
    cache[id] += static_cast<thrust::complex<T>>(A_[(a * nchol + i) * nocc + b] * A_[(b * nchol + i) * nocc + a]);
  }

  __syncthreads(); // required because later on the current thread is accessing
                   // data written by another thread
  int i = nid / 2;
  while (i > 0)
  {
    if (id < i)
      cache[id] += cache[id + i];
    __syncthreads();
    i /= 2; //not sure bitwise operations are actually faster
  }
  if (id == 0)
  {
    T re   = (alp * cache[0]).real();
    T im   = (alp * cache[0]).imag();
    T* re_ = reinterpret_cast<T*>(y + blockIdx.x * incy);
    atomicAdd(re_, re);
    atomicAdd(re_ + 1, im);
  }
}

template<typename T, typename T2>
__global__ void kernel_dot_wanb2(int nwalk,
                                 int nocc,
                                 int nchol,
                                 thrust::complex<T2> const alpha,
                                 thrust::complex<T2> const* Tab,
                                 thrust::complex<T>* y,
                                 int incy)
{
  if (blockIdx.x >= nwalk * nocc * nocc)
    return;
  __shared__ thrust::cuda_cub::core::uninitialized_array<thrust::complex<T>, DOT_BLOCK_SIZE> cache;
  int nocc2              = nocc * nocc;
  int w                  = blockIdx.x / (nocc2);
  int a                  = (blockIdx.x % (nocc2)) / nocc;
  int b                  = (blockIdx.x % (nocc2)) % nocc;
  int i                  = threadIdx.x;
  thrust::complex<T> alp = static_cast<thrust::complex<T>>(alpha);
  thrust::complex<T2> const* A_(Tab + ((w * nocc + a) * nocc) * nchol + b);
  thrust::complex<T2> const* B_(Tab + ((w * nocc + b) * nocc) * nchol + a);
  cache[threadIdx.x] = thrust::complex<T>(0.0);
  while (i < nchol)
  {
    cache[threadIdx.x] += static_cast<thrust::complex<T>>(A_[i * nocc] * B_[i * nocc]);
    i += blockDim.x;
  }
  __syncthreads(); // required because later on the current thread is accessing
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
  {
    T re   = (alp * cache[0]).real();
    T im   = (alp * cache[0]).imag();
    T* re_ = reinterpret_cast<T*>(y + w * incy);
    atomicAdd(re_, re);
    atomicAdd(re_ + 1, im);
  }
}

template<typename T, typename T2>
__global__ void kernel_dot_wpan_waqn_Fwpq(int nwalk,
                                          int nmo,
                                          int nchol,
                                          thrust::complex<T2> const alpha,
                                          thrust::complex<T2> const* Tab,
                                          thrust::complex<T>* F)
{
  __shared__ thrust::cuda_cub::core::uninitialized_array<thrust::complex<T>, DOT_BLOCK_SIZE> cache;
  int p                  = blockIdx.x;
  int q                  = blockIdx.y;
  int a                  = blockIdx.z;
  thrust::complex<T> alp = static_cast<thrust::complex<T>>(alpha);
  for (int w = 0; w < nwalk; w++)
  {
    thrust::complex<T2> const* A_(Tab + ((w * nmo + p) * nmo + a) * nchol);
    thrust::complex<T2> const* B_(Tab + ((w * nmo + a) * nmo + q) * nchol);
    cache[threadIdx.x] = thrust::complex<T>(0.0);
    int i              = threadIdx.x;
    while (i < nchol)
    {
      cache[threadIdx.x] += static_cast<thrust::complex<T>>(A_[i] * B_[i]);
      i += blockDim.x;
    }
    __syncthreads(); // required because later on the current thread is accessing
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
    {
      T re   = (alp * cache[0]).real();
      T im   = (alp * cache[0]).imag();
      T* re_ = reinterpret_cast<T*>(F + (w * nmo + p) * nmo + q);
      atomicAdd(re_, re);
      atomicAdd(re_ + 1, im);
    }
  }
}

void dot_wabn(int nwalk,
              int nocc,
              int nchol,
              std::complex<double> const alpha,
              std::complex<double> const* Tab,
              std::complex<double>* y,
              int incy)
{
  int n_ = nwalk * nocc * nocc;
  dim3 grid_dim(n_, 1, 1);
  kernel_dot_wabn<<<grid_dim, DOT_BLOCK_SIZE>>>(nwalk, nocc, nchol, static_cast<thrust::complex<double> const>(alpha),
                                                reinterpret_cast<thrust::complex<double> const*>(Tab),
                                                reinterpret_cast<thrust::complex<double>*>(y), incy);
  qmc_cuda::cuda_check(cudaGetLastError(), "dot_wabn");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "dot_wabn");
}

void dot_wabn(int nwalk,
              int nocc,
              int nchol,
              std::complex<float> const alpha,
              std::complex<float> const* Tab,
              std::complex<float>* y,
              int incy)
{
  int n_ = nwalk * nocc * nocc;
  dim3 grid_dim(n_, 1, 1);
  kernel_dot_wabn<<<grid_dim, DOT_BLOCK_SIZE>>>(nwalk, nocc, nchol, static_cast<thrust::complex<float> const>(alpha),
                                                reinterpret_cast<thrust::complex<float> const*>(Tab),
                                                reinterpret_cast<thrust::complex<float>*>(y), incy);
  qmc_cuda::cuda_check(cudaGetLastError(), "dot_wabn");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "dot_wabn");
}

void dot_wabn(int nwalk,
              int nocc,
              int nchol,
              std::complex<float> const alpha,
              std::complex<float> const* Tab,
              std::complex<double>* y,
              int incy)
{
  int n_ = nwalk * nocc * nocc;
  dim3 grid_dim(n_, 1, 1);
  kernel_dot_wabn<<<grid_dim, DOT_BLOCK_SIZE>>>(nwalk, nocc, nchol, static_cast<thrust::complex<float> const>(alpha),
                                                reinterpret_cast<thrust::complex<float> const*>(Tab),
                                                reinterpret_cast<thrust::complex<double>*>(y), incy);
  qmc_cuda::cuda_check(cudaGetLastError(), "dot_wabn");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "dot_wabn");
}

// v2
void dot_wanb(int nwalk,
              int nocc,
              int nchol,
              std::complex<double> const alpha,
              std::complex<double> const* Tab,
              std::complex<double>* y,
              int incy)
{
  int a_(8);
  int nf(8);
  int b_ = 1024 / a_;
  int na = (nocc - 1) / a_ + 1;
  int nb = (nocc * nchol - 1) / (b_ * nf) + 1;
  dim3 grid_dim(nwalk, na, nb);
  dim3 block_dim(a_, b_, 1);
  kernel_dot_wanb<<<grid_dim, block_dim>>>(nf, nwalk, nocc, nchol, static_cast<thrust::complex<double> const>(alpha),
                                           reinterpret_cast<thrust::complex<double> const*>(Tab),
                                           reinterpret_cast<thrust::complex<double>*>(y), incy);
  qmc_cuda::cuda_check(cudaGetLastError(), "dot_wanb");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "dot_wanb");
}

void dot_wanb(int nwalk,
              int nocc,
              int nchol,
              std::complex<float> const alpha,
              std::complex<float> const* Tab,
              std::complex<float>* y,
              int incy)
{
  int a_(8);
  int nf(8);
  int b_ = 1024 / a_;
  int na = (nocc - 1) / a_ + 1;
  int nb = (nocc * nchol - 1) / (b_ * nf) + 1;
  dim3 grid_dim(nwalk, na, nb);
  dim3 block_dim(a_, b_, 1);
  kernel_dot_wanb<<<grid_dim, block_dim>>>(nf, nwalk, nocc, nchol, static_cast<thrust::complex<float> const>(alpha),
                                           reinterpret_cast<thrust::complex<float> const*>(Tab),
                                           reinterpret_cast<thrust::complex<float>*>(y), incy);
  qmc_cuda::cuda_check(cudaGetLastError(), "dot_wanb");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "dot_wanb");
}

void dot_wanb(int nwalk,
              int nocc,
              int nchol,
              std::complex<float> const alpha,
              std::complex<float> const* Tab,
              std::complex<double>* y,
              int incy)
{
  int a_(8);
  int nf(8);
  int b_ = 1024 / a_;
  int na = (nocc - 1) / a_ + 1;
  int nb = (nocc * nchol - 1) / (b_ * nf) + 1;
  dim3 grid_dim(nwalk, na, nb);
  dim3 block_dim(a_, b_, 1);
  kernel_dot_wanb<<<grid_dim, block_dim>>>(nf, nwalk, nocc, nchol, static_cast<thrust::complex<float> const>(alpha),
                                           reinterpret_cast<thrust::complex<float> const*>(Tab),
                                           reinterpret_cast<thrust::complex<double>*>(y), incy);
  qmc_cuda::cuda_check(cudaGetLastError(), "dot_wanb");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "dot_wanb");
}

/*
// v2
void dot_wanb( int nwalk, int nocc, int nchol, 
               std::complex<double> const alpha, std::complex<double> const* Tab,
               std::complex<double>* y, int incy)
{
  int n_=nwalk*nocc*nocc;
  dim3 grid_dim(n_,1,1);
  kernel_dot_wanb2<<<grid_dim,DOT_BLOCK_SIZE>>>(nwalk,nocc,nchol,
                                   static_cast<thrust::complex<double> const>(alpha),
                                   reinterpret_cast<thrust::complex<double> const*>(Tab),
                                   reinterpret_cast<thrust::complex<double> *>(y),incy);
  qmc_cuda::cuda_check(cudaGetLastError(),"dot_wanb");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(),"dot_wanb");
}

void dot_wanb( int nwalk, int nocc, int nchol, 
               std::complex<float> const alpha, std::complex<float> const* Tab,
               std::complex<float>* y, int incy)
{
  int n_=nwalk*nocc*nocc;
  dim3 grid_dim(n_,1,1);
  kernel_dot_wanb2<<<grid_dim,DOT_BLOCK_SIZE>>>(nwalk,nocc,nchol,
                                   static_cast<thrust::complex<float> const>(alpha),
                                   reinterpret_cast<thrust::complex<float> const*>(Tab),
                                   reinterpret_cast<thrust::complex<float> *>(y),incy);
  qmc_cuda::cuda_check(cudaGetLastError(),"dot_wanb");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(),"dot_wanb");
}

void dot_wanb( int nwalk, int nocc, int nchol, 
               std::complex<float> const alpha, std::complex<float> const* Tab,
               std::complex<double>* y, int incy)
{
  int n_=nwalk*nocc*nocc;
  dim3 grid_dim(n_,1,1);
  kernel_dot_wanb2<<<grid_dim,DOT_BLOCK_SIZE>>>(nwalk,nocc,nchol,
                                   static_cast<thrust::complex<float> const>(alpha),
                                   reinterpret_cast<thrust::complex<float> const*>(Tab),
                                   reinterpret_cast<thrust::complex<double> *>(y),incy);
  qmc_cuda::cuda_check(cudaGetLastError(),"dot_wanb");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(),"dot_wanb");
}
*/

void dot_wpan_waqn_Fwpq(int nwalk,
                        int nmo,
                        int nchol,
                        std::complex<double> const alpha,
                        std::complex<double> const* Tab,
                        std::complex<double>* F)
{
  dim3 grid_dim(nmo, nmo, nmo);
  kernel_dot_wpan_waqn_Fwpq<<<grid_dim, DOT_BLOCK_SIZE>>>(nwalk, nmo, nchol,
                                                          static_cast<thrust::complex<double> const>(alpha),
                                                          reinterpret_cast<thrust::complex<double> const*>(Tab),
                                                          reinterpret_cast<thrust::complex<double>*>(F));
  qmc_cuda::cuda_check(cudaGetLastError(), "dot_wpan_waqn_Fwpq");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "dot_wpan_waqn_Fwpq");
}

void dot_wpan_waqn_Fwpq(int nwalk,
                        int nmo,
                        int nchol,
                        std::complex<float> const alpha,
                        std::complex<float> const* Tab,
                        std::complex<float>* F)
{
  dim3 grid_dim(nmo, nmo, nmo);
  kernel_dot_wpan_waqn_Fwpq<<<grid_dim, DOT_BLOCK_SIZE>>>(nwalk, nmo, nchol,
                                                          static_cast<thrust::complex<float> const>(alpha),
                                                          reinterpret_cast<thrust::complex<float> const*>(Tab),
                                                          reinterpret_cast<thrust::complex<float>*>(F));
  qmc_cuda::cuda_check(cudaGetLastError(), "dot_wpan_waqn_Fwpq");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "dot_wpan_waqn_Fwpq");
}


void dot_wpan_waqn_Fwpq(int nwalk,
                        int nmo,
                        int nchol,
                        std::complex<float> const alpha,
                        std::complex<float> const* Tab,
                        std::complex<double>* F)
{
  dim3 grid_dim(nmo, nmo, nmo);
  kernel_dot_wpan_waqn_Fwpq<<<grid_dim, DOT_BLOCK_SIZE>>>(nwalk, nmo, nchol,
                                                          static_cast<thrust::complex<float> const>(alpha),
                                                          reinterpret_cast<thrust::complex<float> const*>(Tab),
                                                          reinterpret_cast<thrust::complex<double>*>(F));
  qmc_cuda::cuda_check(cudaGetLastError(), "dot_wpan_waqn_Fwpq");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "dot_wpan_waqn_Fwpq");
}


} // namespace kernels
