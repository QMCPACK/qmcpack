#include "hip/hip_runtime.h"
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

#include<cassert>
#include <complex>
#include<hip/hip_runtime.h>
#include <thrust/complex.h>
#include<hip/hip_runtime.h>
#define ENABLE_HIP 1
#include "AFQMC/Memory/HIP/hip_utilities.h"

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
__global__ void kernel_adotpby(int N, T const alpha, T const* x, int const incx,
                                      T const* y, int const incy, Q const beta, Q* res) {
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
     *res = static_cast<Q>(tmp[0]) + beta *(*res);
   }
   __syncthreads();
}

void adotpby(int N, double const alpha, double const* x, int const incx,
                    double const* y, int const incy,
                    double const beta, double* res)
{
  hipLaunchKernelGGL(kernel_adotpby, dim3(1), dim3(256), 0, 0, N,alpha,x,incx,y,incy,beta,res);
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

void adotpby(int N, std::complex<double> const alpha,
                    std::complex<double> const* x, int const incx,
                    std::complex<double> const* y, int const incy,
                    std::complex<double> const beta, std::complex<double>* res)
{
  hipLaunchKernelGGL(kernel_adotpby, dim3(1), dim3(256), 0, 0, N,
                            static_cast<thrust::complex<double> const >(alpha),
                            reinterpret_cast<thrust::complex<double> const*>(x),incx,
                            reinterpret_cast<thrust::complex<double> const*>(y),incy,
                            static_cast<thrust::complex<double> const>(beta),
                            reinterpret_cast<thrust::complex<double> *>(res));
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

void adotpby(int N, float const alpha, float const* x, int const incx,
                    float const* y, int const incy,
                    float const beta, float* res)
{
  hipLaunchKernelGGL(kernel_adotpby, dim3(1), dim3(256), 0, 0, N,alpha,x,incx,y,incy,beta,res);
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

void adotpby(int N, std::complex<float> const alpha,
                    std::complex<float> const* x, int const incx,
                    std::complex<float> const* y, int const incy,
                    std::complex<float> const beta, std::complex<float>* res)
{
  hipLaunchKernelGGL(kernel_adotpby, dim3(1), dim3(256), 0, 0, N,
                            static_cast<thrust::complex<float> const>(alpha),
                            reinterpret_cast<thrust::complex<float> const*>(x),incx,
                            reinterpret_cast<thrust::complex<float> const*>(y),incy,
                            static_cast<thrust::complex<float> const>(beta),
                            reinterpret_cast<thrust::complex<float> *>(res));
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

void adotpby(int N, float const alpha, float const* x, int const incx,
                    float const* y, int const incy,
                    double const beta, double* res)
{
  hipLaunchKernelGGL(kernel_adotpby, dim3(1), dim3(256), 0, 0, N,alpha,x,incx,y,incy,beta,res);
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

void adotpby(int N, std::complex<float> const alpha,
                    std::complex<float> const* x, int const incx,
                    std::complex<float> const* y, int const incy,
                    std::complex<double> const beta, std::complex<double>* res)
{
  hipLaunchKernelGGL(kernel_adotpby, dim3(1), dim3(256), 0, 0, N,
                            static_cast<thrust::complex<float> const>(alpha),
                            reinterpret_cast<thrust::complex<float> const*>(x),incx,
                            reinterpret_cast<thrust::complex<float> const*>(y),incy,
                            static_cast<thrust::complex<double> const>(beta),
                            reinterpret_cast<thrust::complex<double> *>(res));
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

}

