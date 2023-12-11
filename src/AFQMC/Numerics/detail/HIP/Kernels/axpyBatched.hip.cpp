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
  hipMalloc((void**)&a_, batchSize * sizeof(*a_));
  hipMalloc((void**)&b_, batchSize * sizeof(*b_));
  hipMalloc((void**)&x_, batchSize * sizeof(*x_));
  hipMemcpy(a_, a, batchSize * sizeof(*a), hipMemcpyHostToDevice);
  hipMemcpy(b_, b, batchSize * sizeof(*b), hipMemcpyHostToDevice);
  hipMemcpy(x_, x, batchSize * sizeof(*x), hipMemcpyHostToDevice);
  hipLaunchKernelGGL(kernel_axpy_batched, dim3(batchSize), dim3(128), 0, 0, n, x_, a_, inca, b_, incb, batchSize);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
  hipFree(a_);
  hipFree(b_);
  hipFree(x_);
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
  hipMalloc((void**)&a_, batchSize * sizeof(*a_));
  hipMalloc((void**)&b_, batchSize * sizeof(*b_));
  hipMalloc((void**)&x_, batchSize * sizeof(*x_));
  hipMemcpy(a_, a, batchSize * sizeof(*a), hipMemcpyHostToDevice);
  hipMemcpy(b_, b, batchSize * sizeof(*b), hipMemcpyHostToDevice);
  hipMemcpy(x_, x, batchSize * sizeof(*x), hipMemcpyHostToDevice);
  int nb_(nw > batchSize ? batchSize : nw);
  hipLaunchKernelGGL(kernel_sumGw_batched, dim3(nb_), dim3(256), 0, 0, n, x_, a_, inca, b_, incb, b0, nw, batchSize);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
  hipFree(a_);
  hipFree(b_);
  hipFree(x_);
}

} // namespace kernels
