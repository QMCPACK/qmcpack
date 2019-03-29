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
#include<cuda.h>
#include <thrust/complex.h>
#include<cuda_runtime.h>
#define ENABLE_CUDA 1
#include "AFQMC/Memory/CUDA/cuda_utilities.h"

namespace kernels 
{



template<typename T>
__global__ void kernel_setIdentity(int n, T * A, int lda)
{
  int i = threadIdx.x + blockDim.x*blockIdx.x;
  int j = threadIdx.y + blockDim.y*blockIdx.y;
  if( (i < n) && (j < n) ) 
    if(i==j) {
      A[i*lda+i] = T(1.0); 
    } else {  
      A[j*lda+i] = T(0.0);
    }
}

template<typename T>
__global__ void kernel_setIdentity(int n, thrust::complex<T> * A, int lda)
{
  int i = threadIdx.x + blockDim.x*blockIdx.x;
  int j = threadIdx.y + blockDim.y*blockIdx.y;
  if( (i < n) && (j < n) )
    if(i==j) {
      A[i*lda+i] = thrust::complex<T>(1.0,0.0);
    } else {
      A[j*lda+i] = thrust::complex<T>(0.0,0.0);
    }
}

void setIdentity(int n, double * A, int lda)
{
  int xblock_dim = 16;
  int xgrid_dim = (n + xblock_dim - 1)/xblock_dim;
  dim3 block_dim(xblock_dim,xblock_dim);
  dim3 grid_dim(xgrid_dim,xgrid_dim); 
  kernel_setIdentity<<<grid_dim, block_dim>>>(n,A,lda);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void setIdentity(int n, float * A, int lda)
{
  int xblock_dim = 16;
  int xgrid_dim = (n + xblock_dim - 1)/xblock_dim;
  dim3 block_dim(xblock_dim,xblock_dim);
  dim3 grid_dim(xgrid_dim,xgrid_dim);
  kernel_setIdentity<<<grid_dim, block_dim>>>(n,A,lda);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void setIdentity(int n, std::complex<double> * A, int lda)
{
  int xblock_dim = 16;
  int xgrid_dim = (n + xblock_dim - 1)/xblock_dim;
  dim3 block_dim(xblock_dim,xblock_dim);
  dim3 grid_dim(xgrid_dim,xgrid_dim);
  kernel_setIdentity<<<grid_dim, block_dim>>>(n,reinterpret_cast<thrust::complex<double> *>(A),lda);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void setIdentity(int n, std::complex<float> * A, int lda)
{
  int xblock_dim = 16;
  int xgrid_dim = (n + xblock_dim - 1)/xblock_dim;
  dim3 block_dim(xblock_dim,xblock_dim);
  dim3 grid_dim(xgrid_dim,xgrid_dim);
  kernel_setIdentity<<<grid_dim, block_dim>>>(n,reinterpret_cast<thrust::complex<float> *>(A),lda);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

}

