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
#include "Numerics/detail/cuda_utilities.hpp"

namespace kernels 
{

// Meant to be run with 1 block
template<typename T>
__global__ void kernel_determinant_from_getrf(int N, T *m, int lda, int *piv, T *det) {
   // assert(blockIdx.x==0 and blockIdx.y==0 and blockIdx.z==0)

   __shared__ T tmp[256];
   int t = threadIdx.x;

   tmp[t]=T(1.0);
   int nloop = (N+blockDim.x-1)/blockDim.x;

   for(int i=0, ip=threadIdx.x; i<nloop; i++, ip+=blockDim.x)
    if(ip < N)
    {
      if(piv[ip]==(ip+1)){
        tmp[t] = tmp[t] * m[ip*lda+ip];
      }else{
        tmp[t] = tmp[t] * (-m[ip*lda+ip]);
      }
    }
   __syncthreads();

   // not optimal but ok for now
   if (threadIdx.x == 0) {
     int imax = (N > blockDim.x)?blockDim.x:N;
     for(int i=1; i<imax; i++)
       tmp[0] = tmp[0] * tmp[i];
     *det = tmp[0];
   }
   __syncthreads();
}

void determinant_from_getrf_gpu(int N, double *m, int lda, int *piv, double* res)
{
  kernel_determinant_from_getrf<<<1,256>>>(N,m,lda,piv,res);
  cuda::cuda_check(cudaDeviceSynchronize());
}

void determinant_from_getrf_gpu(int N, std::complex<double> *m, int lda, int *piv, std::complex<double>* res)
{
  kernel_determinant_from_getrf<<<1,256>>>(N,
                                    reinterpret_cast<thrust::complex<double> *>(m),lda,piv,
                                    reinterpret_cast<thrust::complex<double> *>(res) );
  cuda::cuda_check(cudaDeviceSynchronize());
}

}

