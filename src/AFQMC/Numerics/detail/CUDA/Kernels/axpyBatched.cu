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
#define ENABLE_CUDA 1
#include "AFQMC/Memory/CUDA/cuda_utilities.h"

namespace kernels
{
template<typename T>
__global__ void kernel_axpy_batched(int n,
                                    thrust::complex<T>* x,
                                    thrust::complex<T>** a,
                                    int inca,
                                    thrust::complex<T>** b,
                                    int incb,
                                    int batchSize)
{
  int batch = blockIdx.x;
  if (batch >= batchSize)
    return;

  thrust::complex<T>* a_(a[batch]);
  thrust::complex<T>* b_(b[batch]);
  thrust::complex<T> x_(x[batch]);

  int i = threadIdx.x;
  while (i < n)
  {
    b_[i * incb] = b_[i * incb] + x_ * a_[i * inca];
    i += blockDim.x;
  }
}

template<typename T>
__global__ void kernel_sumGw_batched(int n,
                                     thrust::complex<T>* x,
                                     thrust::complex<T>** a,
                                     int inca,
                                     thrust::complex<T>** b,
                                     int incb,
                                     int b0,
                                     int nw,
                                     int batchSize)
{
  if (blockIdx.x >= batchSize)
    return;

  int my_iw = (b0 + blockIdx.x) % nw;

  for (int m = 0; m < batchSize; ++m)
  {
    if ((b0 + m) % nw != my_iw)
      continue;

    thrust::complex<T>* a_(a[m]);
    thrust::complex<T>* b_(b[m]);
    thrust::complex<T> x_(x[m]);

    int i = threadIdx.x;
    while (i < n)
    {
      b_[i * incb] = b_[i * incb] + x_ * a_[i * inca];
      i += blockDim.x;
    }
  }
}

void axpy_batched_gpu(int n,
                      std::complex<double>* x,
                      const std::complex<double>** a,
                      int inca,
                      std::complex<double>** b,
                      int incb,
                      int batchSize)
{
  thrust::complex<double>*x_, **a_, **b_;
  cudaMalloc((void**)&a_, batchSize * sizeof(*a_));
  cudaMalloc((void**)&b_, batchSize * sizeof(*b_));
  cudaMalloc((void**)&x_, batchSize * sizeof(*x_));
  cudaMemcpy(a_, a, batchSize * sizeof(*a), cudaMemcpyHostToDevice);
  cudaMemcpy(b_, b, batchSize * sizeof(*b), cudaMemcpyHostToDevice);
  cudaMemcpy(x_, x, batchSize * sizeof(*x), cudaMemcpyHostToDevice);
  kernel_axpy_batched<<<batchSize, 128>>>(n, x_, a_, inca, b_, incb, batchSize);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
  cudaFree(a_);
  cudaFree(b_);
  cudaFree(x_);
}

void sumGw_batched_gpu(int n,
                       std::complex<double>* x,
                       const std::complex<double>** a,
                       int inca,
                       std::complex<double>** b,
                       int incb,
                       int b0,
                       int nw,
                       int batchSize)
{
  thrust::complex<double>*x_, **a_, **b_;
  cudaMalloc((void**)&a_, batchSize * sizeof(*a_));
  cudaMalloc((void**)&b_, batchSize * sizeof(*b_));
  cudaMalloc((void**)&x_, batchSize * sizeof(*x_));
  cudaMemcpy(a_, a, batchSize * sizeof(*a), cudaMemcpyHostToDevice);
  cudaMemcpy(b_, b, batchSize * sizeof(*b), cudaMemcpyHostToDevice);
  cudaMemcpy(x_, x, batchSize * sizeof(*x), cudaMemcpyHostToDevice);
  int nb_(nw > batchSize ? batchSize : nw);
  kernel_sumGw_batched<<<nb_, 256>>>(n, x_, a_, inca, b_, incb, b0, nw, batchSize);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
  cudaFree(a_);
  cudaFree(b_);
  cudaFree(x_);
}

} // namespace kernels
