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

#ifndef AFQMC_SAMPLEGAUSSIANRNG_H
#define AFQMC_SAMPLEGAUSSIANRNG_H

#include <cassert>
#include <complex>
#include "curand.h"

namespace kernels
{
void sampleGaussianRNG(double* V, int n, curandGenerator_t& gen);
void sampleGaussianRNG(float* V, int n, curandGenerator_t& gen);
void sampleGaussianRNG(std::complex<double>* V, int n, curandGenerator_t& gen);
void sampleGaussianRNG(std::complex<float>* V, int n, curandGenerator_t& gen);

} // namespace kernels

#endif
