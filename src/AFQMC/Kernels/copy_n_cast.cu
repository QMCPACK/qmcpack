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
#define QMC_CUDA 1
#include "AFQMC/Memory/CUDA/cuda_utilities.hpp"

namespace kernels 
{



template<typename T, typename Q>
__global__ void kernel_copy_n_cast(T const* A, int n, Q * B)
{
  int i = threadIdx.x + blockDim.x*blockIdx.x;
  if( i < n ) B[i] = static_cast<Q>(A[i]); 
}

template<typename T, typename Q>
__global__ void kernel_copy_n_cast(thrust::complex<T> const* A, int n, thrust::complex<Q>* B)
{
  int i = threadIdx.x + blockDim.x*blockIdx.x;
  if( i < n ) B[i] = static_cast<thrust::complex<Q> >(A[i]);
}

void copy_n_cast(double const* A, int n, float* B)
{
  int block_dim = 256;
  int grid_dim = (n + block_dim - 1)/block_dim;
  kernel_copy_n_cast<<<grid_dim, block_dim>>>(A,n,B);
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void copy_n_cast(float const* A, int n, double* B)
{
  int block_dim = 256;
  int grid_dim = (n + block_dim - 1)/block_dim;
  kernel_copy_n_cast<<<grid_dim, block_dim>>>(A,n,B);
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void copy_n_cast(std::complex<double> const* A, int n, std::complex<float>* B)
{
  int block_dim = 256;
  int grid_dim = (n + block_dim - 1)/block_dim;
  kernel_copy_n_cast<<<grid_dim, block_dim>>>(reinterpret_cast<thrust::complex<double> const*>(A),n,
                                              reinterpret_cast<thrust::complex<float> *>(B));
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void copy_n_cast(std::complex<float> const* A, int n, std::complex<double>* B)
{
  int block_dim = 256;
  int grid_dim = (n + block_dim - 1)/block_dim;
  kernel_copy_n_cast<<<grid_dim, block_dim>>>(reinterpret_cast<thrust::complex<float> const*>(A),n,
                                              reinterpret_cast<thrust::complex<double> *>(B));
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

}

