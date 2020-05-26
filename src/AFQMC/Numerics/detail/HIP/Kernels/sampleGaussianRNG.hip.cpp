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
#include<hip/hip_runtime.h>
#include "hiprand.h"
#include<hip/hip_runtime.h>
#include "AFQMC/Numerics/detail/HIP/Kernels/hip_settings.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/zero_complex_part.cuh"
#define ENABLE_HIP 1
#include "AFQMC/Memory/HIP/hip_utilities.h"

namespace kernels
{

void sampleGaussianRNG( double* V, int n, hiprandGenerator_t & gen) 
{
  qmc_hip::curand_check(hiprandGenerateNormalDouble(gen,V,n,0.0,1.0),
                                          "hiprandGenerateNormalDouble");
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

 // Convert to double if really necessary
void sampleGaussianRNG( float* V, int n, hiprandGenerator_t & gen) 
{
  qmc_hip::curand_check(hiprandGenerateNormal(gen,V,n,float(0.0),float(1.0)),
                                          "hiprandGenerateNormal");
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

void sampleGaussianRNG( std::complex<double>* V, int n, hiprandGenerator_t & gen) 
{
  qmc_hip::curand_check(hiprandGenerateNormalDouble(gen,
                        reinterpret_cast<double*>(V),2*n,0.0,1.0),
                                          "hiprandGenerateNormalDouble");
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
  // hack hack hack!!!
  kernels::zero_complex_part(n,V);
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
}

void sampleGaussianRNG( std::complex<float>* V, int n, hiprandGenerator_t & gen) 
{
  qmc_hip::curand_check(hiprandGenerateNormal(gen,
                        reinterpret_cast<float*>(V),2*n,float(0.0),float(1.0)),
                                          "hiprandGenerateNormal");
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
  // hack hack hack!!!
  kernels::zero_complex_part(n,V);
  qmc_hip::hip_check(hipGetLastError());
  qmc_hip::hip_check(hipDeviceSynchronize());
} 

}
