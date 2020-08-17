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
#define ENABLE_CUDA 1
#include "AFQMC/Memory/CUDA/cuda_utilities.h"

namespace kernels
{
// Meant to be run with 1 block
/*
template<typename T>
__global__ void kernel_adotpby(int N, T const alpha, T const* x, int const incx, 
                                      T const* y, int const incy, T const beta, T* res) {
   // assert(blockIdx.x==0 and blockIdx.y==0 and blockIdx.z==0)

   __shared__ T tmp[256];
   int t = threadIdx.x;

   tmp[t]=T(0.0);
   int nloop = (N+blockDim.x-1)/blockDim.x;

   for(int i=0, ip=threadIdx.x; i<nloop; i++, ip+=blockDim.x)
    if(ip < N)
    {
      tmp[t] += x[ip*incx]*y[ip*incy]; 
    }
   tmp[t] *= alpha;
   __syncthreads();

   // not optimal but ok for now
   if (threadIdx.x == 0) {
     int imax = (N > blockDim.x)?blockDim.x:N;
     for(int i=1; i<imax; i++)
       tmp[0] += tmp[i];
     *res = tmp[0] + beta *(*res);
   }
   __syncthreads();
}
*/

template<typename T, typename Q>
__global__ void kernel_adotpby(int N,
                               T const alpha,
                               T const* x,
                               int const incx,
                               T const* y,
                               int const incy,
                               Q const beta,
                               Q* res)
{
  // assert(blockIdx.x==0 and blockIdx.y==0 and blockIdx.z==0)

  __shared__ thrust::cuda_cub::core::uninitialized_array<T, 1024> tmp;
  int t = threadIdx.x;

  tmp[t]    = T(0.0);
  int nloop = (N + blockDim.x - 1) / blockDim.x;

  for (int i = 0, ip = threadIdx.x; i < nloop; i++, ip += blockDim.x)
    if (ip < N)
    {
      tmp[t] += x[ip * incx] * y[ip * incy];
    }
  tmp[t] *= alpha;
  __syncthreads();

  t = blockDim.x / 2;
  while (t > 0)
  {
    if (threadIdx.x < t)
      tmp[threadIdx.x] += tmp[threadIdx.x + t];
    __syncthreads();
    t /= 2; //not sure bitwise operations are actually faster
  }
  if (threadIdx.x == 0)
    *res = static_cast<Q>(tmp[0]) + beta * (*res);
}

template<typename T, typename Q>
__global__ void kernel_strided_adotpby(int NB,
                                       int N,
                                       T const alpha,
                                       T const* x,
                                       int const ldx,
                                       T const* y,
                                       int const ldy,
                                       Q const beta,
                                       Q* res,
                                       int inc)
{
  int k = blockIdx.x;
  if (k < NB)
  {
    __shared__ thrust::cuda_cub::core::uninitialized_array<T, 1024> tmp;
    int t = threadIdx.x;

    tmp[t] = T(0.0);

    auto x_(x + k * ldx);
    auto y_(y + k * ldy);
    while (t < N)
    {
      tmp[threadIdx.x] += x_[t] * y_[t];
      t += blockDim.x;
    }
    tmp[threadIdx.x] *= alpha;
    __syncthreads();

    t = blockDim.x / 2;
    while (t > 0)
    {
      if (threadIdx.x < t)
        tmp[threadIdx.x] += tmp[threadIdx.x + t];
      __syncthreads();
      t /= 2; //not sure bitwise operations are actually faster
    }
    if (threadIdx.x == 0)
      *(res + k * inc) = static_cast<Q>(tmp[0]) + beta * (*(res + k * inc));
  }
}

void adotpby(int N,
             double const alpha,
             double const* x,
             int const incx,
             double const* y,
             int const incy,
             double const beta,
             double* res)
{
  kernel_adotpby<<<1, 1024>>>(N, alpha, x, incx, y, incy, beta, res);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void adotpby(int N,
             std::complex<double> const alpha,
             std::complex<double> const* x,
             int const incx,
             std::complex<double> const* y,
             int const incy,
             std::complex<double> const beta,
             std::complex<double>* res)
{
  kernel_adotpby<<<1, 1024>>>(N, static_cast<thrust::complex<double> const>(alpha),
                              reinterpret_cast<thrust::complex<double> const*>(x), incx,
                              reinterpret_cast<thrust::complex<double> const*>(y), incy,
                              static_cast<thrust::complex<double> const>(beta),
                              reinterpret_cast<thrust::complex<double>*>(res));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void adotpby(int N,
             float const alpha,
             float const* x,
             int const incx,
             float const* y,
             int const incy,
             float const beta,
             float* res)
{
  kernel_adotpby<<<1, 1024>>>(N, alpha, x, incx, y, incy, beta, res);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void adotpby(int N,
             std::complex<float> const alpha,
             std::complex<float> const* x,
             int const incx,
             std::complex<float> const* y,
             int const incy,
             std::complex<float> const beta,
             std::complex<float>* res)
{
  kernel_adotpby<<<1, 1024>>>(N, static_cast<thrust::complex<float> const>(alpha),
                              reinterpret_cast<thrust::complex<float> const*>(x), incx,
                              reinterpret_cast<thrust::complex<float> const*>(y), incy,
                              static_cast<thrust::complex<float> const>(beta),
                              reinterpret_cast<thrust::complex<float>*>(res));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void adotpby(int N,
             float const alpha,
             float const* x,
             int const incx,
             float const* y,
             int const incy,
             double const beta,
             double* res)
{
  kernel_adotpby<<<1, 1024>>>(N, alpha, x, incx, y, incy, beta, res);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void adotpby(int N,
             std::complex<float> const alpha,
             std::complex<float> const* x,
             int const incx,
             std::complex<float> const* y,
             int const incy,
             std::complex<double> const beta,
             std::complex<double>* res)
{
  kernel_adotpby<<<1, 1024>>>(N, static_cast<thrust::complex<float> const>(alpha),
                              reinterpret_cast<thrust::complex<float> const*>(x), incx,
                              reinterpret_cast<thrust::complex<float> const*>(y), incy,
                              static_cast<thrust::complex<double> const>(beta),
                              reinterpret_cast<thrust::complex<double>*>(res));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void strided_adotpby(int NB,
                     int N,
                     std::complex<double> const alpha,
                     std::complex<double> const* A,
                     int const lda,
                     std::complex<double> const* B,
                     int const ldb,
                     std::complex<double> const beta,
                     std::complex<double>* C,
                     int ldc)
{
  kernel_strided_adotpby<<<NB, 1024>>>(NB, N, static_cast<thrust::complex<double> const>(alpha),
                                       reinterpret_cast<thrust::complex<double> const*>(A), lda,
                                       reinterpret_cast<thrust::complex<double> const*>(B), ldb,
                                       static_cast<thrust::complex<double> const>(beta),
                                       reinterpret_cast<thrust::complex<double>*>(C), ldc);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void strided_adotpby(int NB,
                     int N,
                     std::complex<float> const alpha,
                     std::complex<float> const* A,
                     int const lda,
                     std::complex<float> const* B,
                     int const ldb,
                     std::complex<double> const beta,
                     std::complex<double>* C,
                     int ldc)
{
  kernel_strided_adotpby<<<NB, 1024>>>(NB, N, static_cast<thrust::complex<float> const>(alpha),
                                       reinterpret_cast<thrust::complex<float> const*>(A), lda,
                                       reinterpret_cast<thrust::complex<float> const*>(B), ldb,
                                       static_cast<thrust::complex<double> const>(beta),
                                       reinterpret_cast<thrust::complex<double>*>(C), ldc);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

} // namespace kernels
