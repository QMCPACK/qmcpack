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
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/uninitialized_fill.h>
#define QMC_CUDA 1
#include "Numerics/detail/cuda_utilities.hpp"
//#include "AFQMC/Kernels/strided_range.hpp"

namespace kernels 
{

template<typename T>
__global__ void kernel_uninitialized_fill_n(int N, T* x, T const a) 
{
   int nloop = (N+blockDim.x-1)/blockDim.x;
   for(int i=0, ip=threadIdx.x; i<nloop; i++, ip+=blockDim.x)
    if(ip < N)
    {
      x[ip] = a; 
    }
   __syncthreads();
}

void uninitialized_fill_n(bool * first, int N, bool const value)
{
  kernel_uninitialized_fill_n<<<1,256>>>(N,first,value);
  cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_fill_n(int * first, int N, int const value)
{ 
  kernel_uninitialized_fill_n<<<1,256>>>(N,first,value);
  cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_fill_n(float * first, int N, float const value)
{
  kernel_uninitialized_fill_n<<<1,256>>>(N,first,value);
  cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_fill_n(double * first, int N, double const value)
{
  kernel_uninitialized_fill_n<<<1,256>>>(N,first,value);
  cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_fill_n(std::complex<float> * first, int N, std::complex<float> const value)
{
  kernel_uninitialized_fill_n<<<1,256>>>(N,first,value);
  cuda::cuda_check(cudaDeviceSynchronize());
}

void uninitialized_fill_n(std::complex<double> * first, int N, std::complex<double> const value)
{ 
  kernel_uninitialized_fill_n<<<1,256>>>(N,first,value);
  cuda::cuda_check(cudaDeviceSynchronize());
} 

void uninitialized_fill_n(double2 * first, int N, double2 const value)
{ 
  kernel_uninitialized_fill_n<<<1,256>>>(N,first,value);
  cuda::cuda_check(cudaDeviceSynchronize());
}

}
