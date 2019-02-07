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
#include <type_traits>
#include <thrust/complex.h>
#include <thrust/device_ptr.h>
#include <thrust/fill.h>
#define QMC_CUDA 1
#include "AFQMC/Memory/CUDA/cuda_utilities.h"
//#include "AFQMC/Kernels/strided_range.hpp"

namespace kernels 
{

template<typename T, typename Size>
__global__ void kernel_fill_n(Size N, T* x, int incx, T const a)
{
 for(Size ip=Size(threadIdx.x); ip<N; ip+=Size(blockDim.x))
  {
    x[ip*incx] = a;
  }
}


void fill_n(int * first, int N, int incx, int const value)
{ 
  kernel_fill_n<<<1,256>>>(N,first,1,value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void fill_n(float * first, int N, int incx, float const value)
{ 
  kernel_fill_n<<<1,256>>>(N,first,1,value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void fill_n(double * first, int N, int incx, double const value)
{ 
  kernel_fill_n<<<1,256>>>(N,first,1,value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void fill_n(std::complex<float> * first, int N, int incx, std::complex<float> const value)
{ 
  kernel_fill_n<<<1,256>>>(N,first,1,value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void fill_n(std::complex<double> * first, int N, int stride, std::complex<double> const value)
{ 
  kernel_fill_n<<<1,256>>>(N,first,1,value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}


void fill_n(int * first, int N, int const value)
{ 
  kernel_fill_n<<<1,256>>>(N,first,1,value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void fill_n(float * first, int N, float const value)
{ 
  kernel_fill_n<<<1,256>>>(N,first,1,value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void fill_n(double * first, int N, double const value)
{ 
  kernel_fill_n<<<1,256>>>(N,first,1,value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void fill_n(std::complex<float> * first, int N, std::complex<float> const value)
{ 
  kernel_fill_n<<<1,256>>>(N,first,1,value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void fill_n(std::complex<double> * first, int N, std::complex<double> const value)
{ 
  kernel_fill_n<<<1,256>>>(N,first,1,value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}


}

