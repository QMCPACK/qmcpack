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

#ifndef AFQMC_COPY_N_CAST_KERNELS_HPP
#define AFQMC_COPY_N_CAST_KERNELS_HPP

#include<cassert>
#include <complex>

namespace kernels
{

void copy_n_cast(double const* A, int n, float* B);
void copy_n_cast(float const* A, int n, double* B);
void copy_n_cast(std::complex<double> const* A, int n, std::complex<float>* B);
void copy_n_cast(std::complex<float> const* A, int n, std::complex<double>* B);

}

#endif
