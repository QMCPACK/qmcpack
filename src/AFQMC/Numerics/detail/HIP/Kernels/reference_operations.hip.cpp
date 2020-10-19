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


#include <complex>
#include <thrust/complex.h>
#include "AFQMC/Numerics/detail/HIP/hip_kernel_utils.h"

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
  hipLaunchKernelGGL(op_plus__, dim3(1), dim3(1), 0, 0, x, inc);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void op_plus(float* x, float inc)
{
  hipLaunchKernelGGL(op_plus__, dim3(1), dim3(1), 0, 0, x, inc);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void op_plus(std::complex<double>* x, std::complex<double> inc)
{
  hipLaunchKernelGGL(op_plus__, dim3(1), dim3(1), 0, 0, reinterpret_cast<thrust::complex<double>*>(x),
                     static_cast<thrust::complex<double>>(inc));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void op_plus(std::complex<float>* x, std::complex<float> inc)
{
  hipLaunchKernelGGL(op_plus__, dim3(1), dim3(1), 0, 0, reinterpret_cast<thrust::complex<float>*>(x),
                     static_cast<thrust::complex<float>>(inc));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

// -=
void op_minus(double* x, double inc)
{
  hipLaunchKernelGGL(op_minus__, dim3(1), dim3(1), 0, 0, x, inc);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void op_minus(float* x, float inc)
{
  hipLaunchKernelGGL(op_minus__, dim3(1), dim3(1), 0, 0, x, inc);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void op_minus(std::complex<double>* x, std::complex<double> inc)
{
  hipLaunchKernelGGL(op_minus__, dim3(1), dim3(1), 0, 0, reinterpret_cast<thrust::complex<double>*>(x),
                     static_cast<thrust::complex<double>>(inc));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void op_minus(std::complex<float>* x, std::complex<float> inc)
{
  hipLaunchKernelGGL(op_minus__, dim3(1), dim3(1), 0, 0, reinterpret_cast<thrust::complex<float>*>(x),
                     static_cast<thrust::complex<float>>(inc));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

// *=
void op_times(double* x, double inc)
{
  hipLaunchKernelGGL(op_times__, dim3(1), dim3(1), 0, 0, x, inc);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void op_times(float* x, float inc)
{
  hipLaunchKernelGGL(op_times__, dim3(1), dim3(1), 0, 0, x, inc);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void op_times(std::complex<double>* x, std::complex<double> inc)
{
  hipLaunchKernelGGL(op_times__, dim3(1), dim3(1), 0, 0, reinterpret_cast<thrust::complex<double>*>(x),
                     static_cast<thrust::complex<double>>(inc));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void op_times(std::complex<float>* x, std::complex<float> inc)
{
  hipLaunchKernelGGL(op_times__, dim3(1), dim3(1), 0, 0, reinterpret_cast<thrust::complex<float>*>(x),
                     static_cast<thrust::complex<float>>(inc));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

// /=
void op_div(double* x, double inc)
{
  hipLaunchKernelGGL(op_div__, dim3(1), dim3(1), 0, 0, x, inc);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void op_div(float* x, float inc)
{
  hipLaunchKernelGGL(op_div__, dim3(1), dim3(1), 0, 0, x, inc);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void op_div(std::complex<double>* x, std::complex<double> inc)
{
  hipLaunchKernelGGL(op_div__, dim3(1), dim3(1), 0, 0, reinterpret_cast<thrust::complex<double>*>(x),
                     static_cast<thrust::complex<double>>(inc));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void op_div(std::complex<float>* x, std::complex<float> inc)
{
  hipLaunchKernelGGL(op_div__, dim3(1), dim3(1), 0, 0, reinterpret_cast<thrust::complex<float>*>(x),
                     static_cast<thrust::complex<float>>(inc));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}


} // namespace kernels
