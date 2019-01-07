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
#include "AFQMC/Kernels/cuda_settings.h"
#define QMC_CUDA 1
#include "Numerics/detail/cuda_utilities.hpp"

namespace kernels 
{

template<typename T, typename T2>
__global__ void kernel_dot_wabn_wban( int nw, int na, int nb, int nc, 
                    T const alpha, T2 const* A, T2 const* B, T* y, int incy) 
{
    if( blockIdx.x > nw*na*nb ) return;
    __shared__ T cache[ DOT_BLOCK_SIZE ];
    int w = blockIdx.x/(na*nb);
    int a = (blockIdx.x%(na*nb))/nb;
    int b = (blockIdx.x%(na*nb))%nb;
    int i = threadIdx.x;
    T2 const* A_(A + (w*na+a)*nb + b);
    T2 const* B_(B + (w*nb+b)*na + a);
    cache[ threadIdx.x ] = T(0.0);
    while( i < nc ) {
        cache[ threadIdx.x ] += static_cast<T>(A_[ i ] * B_[ i ]);
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
    //if( threadIdx.x == 0 ) *(y+w*incy) = (*(y+w*incy)) + alpha * cache[ 0 ];
    if( threadIdx.x == 0 ) atomicAdd_system(y+w*incy,alpha * cache[ 0 ]); 
}

template<typename T, typename T2>
__global__ void kernel_dot_wabn_wban( int nw, int na, int nb, int nc,
                    thrust::complex<T> const alpha, thrust::complex<T2> const* A, 
                    thrust::complex<T2> const* B, thrust::complex<T>* y, int incy)
{
    if( blockIdx.x > nw*na*nb ) return;
    __shared__ thrust::complex<T> cache[ DOT_BLOCK_SIZE ];
    int w = blockIdx.x/(na*nb);
    int a = (blockIdx.x%(na*nb))/nb;
    int b = (blockIdx.x%(na*nb))%nb;
    int i = threadIdx.x;
    thrust::complex<T2> const* A_(A + ((w*na+a)*nb + b)*nc);
    thrust::complex<T2> const* B_(B + ((w*nb+b)*na + a)*nc);
    cache[ threadIdx.x ] = thrust::complex<T>(0.0);
    while( i < nc ) {
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
    //if( threadIdx.x == 0 ) *(y+w*incy) = (*(y+w*incy)) + alpha * cache[ 0 ];
    if( threadIdx.x == 0 ) {
        T2 re = (alpha * cache[ 0 ]).real();
        T2 im = (alpha * cache[ 0 ]).imag();
        T2* re_ = reinterpret_cast<T2*>(y+w*incy);
        atomicAdd_system(re_,re); 
        atomicAdd_system(re_+1,im); 
    }
}

void dot_wabn_wban( int nw, int na, int nb, int nc,
                    double const alpha, double const* A, double const* B, double* y, int incy)
{
  int n_=nw*na*nb;
  kernel_dot_wabn_wban<<<n_,DOT_BLOCK_SIZE>>>(nw,na,nb,nc,alpha,A,B,y,incy);
  cuda::cuda_check(cudaDeviceSynchronize());
}

void dot_wabn_wban( int nw, int na, int nb, int nc,
                    double const alpha, float const* A, float const* B, double* y, int incy)
{
  int n_=nw*na*nb;
  kernel_dot_wabn_wban<<<n_,DOT_BLOCK_SIZE>>>(nw,na,nb,nc,alpha,A,B,y,incy);
  cuda::cuda_check(cudaDeviceSynchronize());
}

void dot_wabn_wban( int nw, int na, int nb, int nc,
                    std::complex<double> const alpha, std::complex<double> const* A, 
                    std::complex<double> const* B, std::complex<double>* y, int incy)
{
  int n_=nw*na*nb;
  kernel_dot_wabn_wban<<<n_,DOT_BLOCK_SIZE>>>(nw,na,nb,nc,
                                   static_cast<thrust::complex<double> const>(alpha),
                                   reinterpret_cast<thrust::complex<double> const*>(A),
                                   reinterpret_cast<thrust::complex<double> const*>(B),
                                   reinterpret_cast<thrust::complex<double> *>(y),incy);
  cuda::cuda_check(cudaDeviceSynchronize());
}

void dot_wabn_wban( int nw, int na, int nb, int nc,
                    std::complex<double> const alpha, std::complex<float> const* A, 
                    std::complex<float> const* B, std::complex<double>* y, int incy)
{
  int n_=nw*na*nb;
  kernel_dot_wabn_wban<<<n_,DOT_BLOCK_SIZE>>>(nw,na,nb,nc,
                                   static_cast<thrust::complex<double> const>(alpha),
                                   reinterpret_cast<thrust::complex<float> const*>(A),
                                   reinterpret_cast<thrust::complex<float> const*>(B),
                                   reinterpret_cast<thrust::complex<double> *>(y),incy);
  cuda::cuda_check(cudaDeviceSynchronize());
}

}
