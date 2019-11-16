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
#include "AFQMC/Numerics/detail/CUDA/Kernels/cuda_settings.h"
#define ENABLE_CUDA 1
#include "AFQMC/Memory/CUDA/cuda_utilities.h"
#if __CUDA_ARCH__ < 600
#include "AFQMC/Numerics/detail/CUDA/Kernels/myAtomicAdd.cu"
#endif

namespace kernels 
{

// Tab nwalk][nocc][nocc][nchol]
template<typename T, typename T2>
__global__ void kernel_dot_wabn(int nwalk, int nocc, int nchol,
                    thrust::complex<T2> const alpha, thrust::complex<T2> const* Tab, 
                    thrust::complex<T>* y, int incy)
{
    if( blockIdx.x >= nwalk*nocc*nocc ) return;
    __shared__ thrust::complex<T> cache[ DOT_BLOCK_SIZE ];
    int nocc2 = nocc*nocc;
    int w = blockIdx.x/(nocc2);
    int a = (blockIdx.x%(nocc2))/nocc;
    int b = (blockIdx.x%(nocc2))%nocc;
    int i = threadIdx.x;
    thrust::complex<T> alp = static_cast<thrust::complex<T>>(alpha);
    thrust::complex<T2> const* A_(Tab + ((w*nocc+a)*nocc + b)*nchol);
    thrust::complex<T2> const* B_(Tab + ((w*nocc+b)*nocc + a)*nchol);
    cache[ threadIdx.x ] = thrust::complex<T>(0.0);
    while( i < nchol ) {
        cache[ threadIdx.x ] += static_cast<thrust::complex<T>>(A_[ i ] * B_[ i ]);
        i += blockDim.x;
    }
    __syncthreads(); // required because later on the current thread is accessing
                     // data written by another thread    
    i = DOT_BLOCK_SIZE / 2;
    while( i > 0 ) {
        if( threadIdx.x < i ) cache[ threadIdx.x ] += cache[ threadIdx.x + i ];
        __syncthreads();
        i /= 2; //not sure bitwise operations are actually faster
    }
    if( threadIdx.x == 0 ) {
        T re = (alp * cache[ 0 ]).real();
        T im = (alp * cache[ 0 ]).imag();
        T* re_ = reinterpret_cast<T*>(y+w*incy);
#if __CUDA_ARCH__ < 600
        myAtomicAdd(re_,re); 
        myAtomicAdd(re_+1,im); 
#else
        atomicAdd(re_,re); 
        atomicAdd(re_+1,im); 
#endif
    }
}

void dot_wabn( int nwalk, int nocc, int nchol, 
               std::complex<double> const alpha, std::complex<double> const* Tab, 
               std::complex<double>* y, int incy)
{
  int n_=nwalk*nocc*nocc;
  dim3 grid_dim(n_,1,1);
  kernel_dot_wabn<<<grid_dim,DOT_BLOCK_SIZE>>>(nwalk,nocc,nchol,
                                   static_cast<thrust::complex<double> const>(alpha),
                                   reinterpret_cast<thrust::complex<double> const*>(Tab),
                                   reinterpret_cast<thrust::complex<double> *>(y),incy);
  qmc_cuda::cuda_check(cudaGetLastError(),"dot_wabn");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(),"dot_wabn");
}

void dot_wabn( int nwalk, int nocc, int nchol,                         
               std::complex<float> const alpha, std::complex<float> const* Tab,      
               std::complex<float>* y, int incy)
{
  int n_=nwalk*nocc*nocc;
  dim3 grid_dim(n_,1,1);
  kernel_dot_wabn<<<grid_dim,DOT_BLOCK_SIZE>>>(nwalk,nocc,nchol,
                                   static_cast<thrust::complex<float> const>(alpha),
                                   reinterpret_cast<thrust::complex<float> const*>(Tab),
                                   reinterpret_cast<thrust::complex<float> *>(y),incy);
  qmc_cuda::cuda_check(cudaGetLastError(),"dot_wabn");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(),"dot_wabn");
}

void dot_wabn( int nwalk, int nocc, int nchol,                         
               std::complex<float> const alpha, std::complex<float> const* Tab,      
               std::complex<double>* y, int incy)
{
  int n_=nwalk*nocc*nocc;
  dim3 grid_dim(n_,1,1);
  kernel_dot_wabn<<<grid_dim,DOT_BLOCK_SIZE>>>(nwalk,nocc,nchol,
                                   static_cast<thrust::complex<float> const>(alpha),
                                   reinterpret_cast<thrust::complex<float> const*>(Tab),
                                   reinterpret_cast<thrust::complex<double> *>(y),incy);
  qmc_cuda::cuda_check(cudaGetLastError(),"dot_wabn");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(),"dot_wabn");
}

}
