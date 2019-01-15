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
#include "curand.h"
#include<cuda_runtime.h>
#include "AFQMC/Kernels/cuda_settings.h"
#include "AFQMC/Kernels/zero_complex_part.hpp"
#define QMC_CUDA 1

namespace kernels
{

void sampleGaussianRNG( double* V, int n, curandGenerator_t & gen) 
{
  qmc_cuda::cuda_check(curandGenerateNormalDouble(gen,V,n,0.0,1.0),
                                          "curandGenerateNormalDouble");
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void sampleGaussianRNG( float* V, int n, curandGenerator_t & gen) 
{
  qmc_cuda::cuda_check(curandGenerateNormalFloat(gen,V,n,0.0,1.0),
                                          "curandGenerateNormalFloat");
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void sampleGaussianRNG( std::complex<double>* V, int n, curandGenerator_t & gen) 
{
  qmc_cuda::cuda_check(curandGenerateNormalDouble(gen,
                        reinterpret_cast<double*>(V),2*n,0.0,1.0),
                                          "curandGenerateNormalDouble");
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
  // hack hack hack!!!
  kernels::zero_complex_part(n,V);
}

void sampleGaussianRNG( std::complex<float>* V, int n, curandGenerator_t & gen) 
{
  qmc_cuda::cuda_check(curandGenerateNormalFloat(gen,
                        reinterpret_cast<float*>(V),2*n,0.0,1.0),
                                          "curandGenerateNormalFloat");
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
  // hack hack hack!!!
  kernels::zero_complex_part(n,V);
} 

}
