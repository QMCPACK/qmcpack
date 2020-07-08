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


#include<cassert>
#include <complex>
#include <type_traits>
/*
#include <thrust/complex.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/uninitialized_fill.h>
*/
#include "AFQMC/Memory/HIP/hip_utilities.h"

namespace kernels
{

template<typename T, typename Size>
__global__ void kernel_uninitialized_fill_n(Size N, T* x, T const a)
{
   Size nloop = Size((N+blockDim.x-1)/blockDim.x);
   for(Size i=0, ip=Size(threadIdx.x); i<nloop; i++, ip+=Size(blockDim.x))
    if(ip < N)
    {
      x[ip] = a;
    }
   __syncthreads();
}

void uninitialized_fill_n(bool * first, int N, bool const value)
{
  hipLaunchKernelGGL(kernel_uninitialized_fill_n, dim3(1), dim3(256), 0, 0, N,first,value);
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

void uninitialized_fill_n(int * first, int N, int const value)
{
  hipLaunchKernelGGL(kernel_uninitialized_fill_n, dim3(1), dim3(256), 0, 0, N,first,value);
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

void uninitialized_fill_n(float * first, int N, float const value)
{
  hipLaunchKernelGGL(kernel_uninitialized_fill_n, dim3(1), dim3(256), 0, 0, N,first,value);
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

void uninitialized_fill_n(double * first, int N, double const value)
{
  hipLaunchKernelGGL(kernel_uninitialized_fill_n, dim3(1), dim3(256), 0, 0, N,first,value);
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

void uninitialized_fill_n(std::complex<float> * first, int N, std::complex<float> const value)
{
  hipLaunchKernelGGL(kernel_uninitialized_fill_n, dim3(1), dim3(256), 0, 0, N,first,value);
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

void uninitialized_fill_n(std::complex<double> * first, int N, std::complex<double> const value)
{
  hipLaunchKernelGGL(kernel_uninitialized_fill_n, dim3(1), dim3(256), 0, 0, N,first,value);
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

void uninitialized_fill_n(double2 * first, int N, double2 const value)
{
  hipLaunchKernelGGL(kernel_uninitialized_fill_n, dim3(1), dim3(256), 0, 0, N,first,value);
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

void uninitialized_fill_n(bool * first, long N, bool const value)
{
  hipLaunchKernelGGL(kernel_uninitialized_fill_n, dim3(1), dim3(256), 0, 0, N,first,value);
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

void uninitialized_fill_n(int * first, long N, int const value)
{
  hipLaunchKernelGGL(kernel_uninitialized_fill_n, dim3(1), dim3(256), 0, 0, N,first,value);
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

void uninitialized_fill_n(float * first, long N, float const value)
{
  hipLaunchKernelGGL(kernel_uninitialized_fill_n, dim3(1), dim3(256), 0, 0, N,first,value);
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

void uninitialized_fill_n(double * first, long N, double const value)
{
  hipLaunchKernelGGL(kernel_uninitialized_fill_n, dim3(1), dim3(256), 0, 0, N,first,value);
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

void uninitialized_fill_n(std::complex<float> * first, long N, std::complex<float> const value)
{
  hipLaunchKernelGGL(kernel_uninitialized_fill_n, dim3(1), dim3(256), 0, 0, N,first,value);
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

void uninitialized_fill_n(std::complex<double> * first, long N, std::complex<double> const value)
{
  hipLaunchKernelGGL(kernel_uninitialized_fill_n, dim3(1), dim3(256), 0, 0, N,first,value);
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

void uninitialized_fill_n(double2 * first, long N, double2 const value)
{
  hipLaunchKernelGGL(kernel_uninitialized_fill_n, dim3(1), dim3(256), 0, 0, N,first,value);
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}


}
