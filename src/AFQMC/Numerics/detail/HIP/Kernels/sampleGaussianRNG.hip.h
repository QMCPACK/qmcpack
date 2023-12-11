///////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Fionn Malone, malone14@llnl.gov, LLNL
//
// File created by: Fionn Malone, malone14@llnl.gov, LLNL
////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_SAMPLEGAUSSIANRNG_H
#define AFQMC_SAMPLEGAUSSIANRNG_H

#include <cassert>
#include <complex>
#include <rocrand/rocrand.h>

namespace kernels
{
void sampleGaussianRNG(double* V, int n, rocrand_generator& gen);
void sampleGaussianRNG(float* V, int n, rocrand_generator& gen);
void sampleGaussianRNG(std::complex<double>* V, int n, rocrand_generator& gen);
void sampleGaussianRNG(std::complex<float>* V, int n, rocrand_generator& gen);

} // namespace kernels

#endif
