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
#define QMC_CUDA 1
#include "AFQMC/Memory/CUDA/cuda_utilities.h"
//#include "AFQMC/Kernels/strided_range.hpp"

namespace kernels 
{

template<typename T, typename Size>
__global__ void kernel_uninitialized_copy_n(Size N, T* x, T const* arr)
{
   Size nloop = Size((N+blockDim.x-1)/blockDim.x);
   for(Size i=0, ip=Size(threadIdx.x); i<nloop; i++, ip+=Size(blockDim.x))
    if(ip < N)
    {
      x[ip] = arr[ip];
    }
   __syncthreads();
}


void uninitialized_copy_n(double * first, int N, double const* array)
{
  kernel_uninitialized_copy_n<<<1,256>>>(N,first,array);
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_copy_n(std::complex<double>* first, int N, std::complex<double> const* array)
{
  kernel_uninitialized_copy_n<<<1,256>>>(N,first,array);
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_copy_n(int* first, int N, int const* array)
{
  kernel_uninitialized_copy_n<<<1,256>>>(N,first,array);
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

// long
void uninitialized_copy_n(double * first, long N, double const* array)
{
  kernel_uninitialized_copy_n<<<1,256>>>(N,first,array);
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_copy_n(std::complex<double>* first, long N, std::complex<double> const* array)
{
  kernel_uninitialized_copy_n<<<1,256>>>(N,first,array);
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_copy_n(int* first, long N, int const* array)
{
  kernel_uninitialized_copy_n<<<1,256>>>(N,first,array);
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

}
