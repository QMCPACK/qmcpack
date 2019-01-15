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

#ifndef AFQMC_UNINITIALIZED_COPY_N_KERNELS_HPP
#define AFQMC_UNINITIALIZED_COPY_N_KERNELS_HPP

#include<cassert>
#include <complex>

namespace kernels 
{

void uninitialized_copy_n(double * first, int N, double const* array);
void uninitialized_copy_n(std::complex<double>* first, int N, std::complex<double> const* array);
void uninitialized_copy_n(int* first, int N, int const* array);

// long
void uninitialized_copy_n(double * first, long N, double const* array);
void uninitialized_copy_n(std::complex<double>* first, long N, std::complex<double> const* array);
void uninitialized_copy_n(int* first, long N, int const* array);

}

#endif
