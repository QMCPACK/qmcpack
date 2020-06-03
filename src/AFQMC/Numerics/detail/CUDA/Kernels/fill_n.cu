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

#define BOOST_NO_AUTO_PTR
#include<cassert>
#include <complex>
#include <type_traits>
#define ENABLE_CUDA 1
#include "AFQMC/Memory/CUDA/cuda_utilities.h"

namespace kernels 
{

template<typename Size, typename Size1, typename T>
__global__ void kernel_fill_n(Size N, T* x, Size1 incx, T const a)
{
 if(threadIdx.x >= N) return; 
 for(Size ip=Size(threadIdx.x); ip<N; ip+=Size(blockDim.x))
  {
    x[ip*Size(incx)] = a;
  }
}

template<typename T, typename Size>
__global__ void kernel_fill2D_n(Size N, Size M, T* y, Size lda, T const a)
{
  for(Size ip=Size(threadIdx.x); ip<N; ip+=Size(blockDim.x))
    for(Size jp=Size(threadIdx.y); jp<M; jp+=Size(blockDim.y))
    {
      y[ip*lda + jp] = a;
    }
}

void fill_n(char * first, int N, int incx, char const value)
{
  kernel_fill_n<<<1,256>>>(N,first,incx,value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void fill_n(int * first, int N, int incx, int const value)
{ 
  kernel_fill_n<<<1,256>>>(N,first,incx,value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void fill_n(float * first, int N, int incx, float const value)
{ 
  kernel_fill_n<<<1,256>>>(N,first,incx,value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void fill_n(double * first, int N, int incx, double const value)
{ 
  kernel_fill_n<<<1,256>>>(N,first,incx,value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void fill_n(std::complex<float> * first, int N, int incx, std::complex<float> const value)
{ 
  kernel_fill_n<<<1,256>>>(N,first,incx,value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void fill_n(std::complex<double> * first, int N, int incx, std::complex<double> const value)
{ 
  kernel_fill_n<<<1,256>>>(N,first,incx,value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void fill_n(char * first, int N, char const value)
{
  kernel_fill_n<<<1,256>>>(N,first,1,value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void fill_n(long int * first, long unsigned int N, const long int value)
{
  kernel_fill_n<<<1,256>>>(N,first,1,value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void fill_n(long unsigned int * first, long unsigned int N, const long unsigned int value)
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

void fill2D_n(int N, int M, int* A, int lda, int const value)
{
  kernel_fill2D_n<<<32,32>>>(N,M,A,lda,value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void fill2D_n(int N, int M, float* A, int lda, float const value)
{
  kernel_fill2D_n<<<32,32>>>(N,M,A,lda,value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void fill2D_n(int N, int M, double* A, int lda, double const value)
{
  kernel_fill2D_n<<<32,32>>>(N,M,A,lda,value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void fill2D_n(int N, int M, std::complex<double>* A, int lda, std::complex<double> const value)
{
  kernel_fill2D_n<<<32,32>>>(N,M,A,lda,value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void fill2D_n(int N, int M, std::complex<float>* A, int lda, std::complex<float> const value)
{
  kernel_fill2D_n<<<32,32>>>(N,M,A,lda,value);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

}

