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

#include <complex>
#include <thrust/complex.h>
#define ENABLE_CUDA 1
#include "AFQMC/Memory/CUDA/cuda_utilities.h"

namespace kernels
{
template<typename T>
__global__ void op_plus__(T* x, T inc)
{
  if (threadIdx.x == 0)
    *x += inc;
}

template<typename T>
__global__ void op_plus__(thrust::complex<T>* x, thrust::complex<T> inc)
{
  if (threadIdx.x == 0)
    *x += inc;
}

template<typename T>
__global__ void op_minus__(T* x, T inc)
{
  if (threadIdx.x == 0)
    *x -= inc;
}

template<typename T>
__global__ void op_minus__(thrust::complex<T>* x, thrust::complex<T> inc)
{
  if (threadIdx.x == 0)
    *x -= inc;
}


template<typename T>
__global__ void op_times__(T* x, T inc)
{
  if (threadIdx.x == 0)
    *x *= inc;
}

template<typename T>
__global__ void op_times__(thrust::complex<T>* x, thrust::complex<T> inc)
{
  if (threadIdx.x == 0)
    *x *= inc;
}


template<typename T>
__global__ void op_div__(T* x, T inc)
{
  if (threadIdx.x == 0)
    *x /= inc;
}

template<typename T>
__global__ void op_div__(thrust::complex<T>* x, thrust::complex<T> inc)
{
  if (threadIdx.x == 0)
    *x /= inc;
}


// +=
void op_plus(double* x, double inc)
{
  op_plus__<<<1, 1>>>(x, inc);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void op_plus(float* x, float inc)
{
  op_plus__<<<1, 1>>>(x, inc);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void op_plus(std::complex<double>* x, std::complex<double> inc)
{
  op_plus__<<<1, 1>>>(reinterpret_cast<thrust::complex<double>*>(x), static_cast<thrust::complex<double>>(inc));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void op_plus(std::complex<float>* x, std::complex<float> inc)
{
  op_plus__<<<1, 1>>>(reinterpret_cast<thrust::complex<float>*>(x), static_cast<thrust::complex<float>>(inc));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

// -=
void op_minus(double* x, double inc)
{
  op_minus__<<<1, 1>>>(x, inc);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void op_minus(float* x, float inc)
{
  op_minus__<<<1, 1>>>(x, inc);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void op_minus(std::complex<double>* x, std::complex<double> inc)
{
  op_minus__<<<1, 1>>>(reinterpret_cast<thrust::complex<double>*>(x), static_cast<thrust::complex<double>>(inc));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void op_minus(std::complex<float>* x, std::complex<float> inc)
{
  op_minus__<<<1, 1>>>(reinterpret_cast<thrust::complex<float>*>(x), static_cast<thrust::complex<float>>(inc));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

// *=
void op_times(double* x, double inc)
{
  op_times__<<<1, 1>>>(x, inc);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void op_times(float* x, float inc)
{
  op_times__<<<1, 1>>>(x, inc);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void op_times(std::complex<double>* x, std::complex<double> inc)
{
  op_times__<<<1, 1>>>(reinterpret_cast<thrust::complex<double>*>(x), static_cast<thrust::complex<double>>(inc));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void op_times(std::complex<float>* x, std::complex<float> inc)
{
  op_times__<<<1, 1>>>(reinterpret_cast<thrust::complex<float>*>(x), static_cast<thrust::complex<float>>(inc));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

// /=
void op_div(double* x, double inc)
{
  op_div__<<<1, 1>>>(x, inc);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void op_div(float* x, float inc)
{
  op_div__<<<1, 1>>>(x, inc);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void op_div(std::complex<double>* x, std::complex<double> inc)
{
  op_div__<<<1, 1>>>(reinterpret_cast<thrust::complex<double>*>(x), static_cast<thrust::complex<double>>(inc));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void op_div(std::complex<float>* x, std::complex<float> inc)
{
  op_div__<<<1, 1>>>(reinterpret_cast<thrust::complex<float>*>(x), static_cast<thrust::complex<float>>(inc));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}


} // namespace kernels
