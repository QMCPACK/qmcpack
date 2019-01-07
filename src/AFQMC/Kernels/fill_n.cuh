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

#ifndef AFQMC_FILL_N_KERNELS_HPP
#define AFQMC_FILL_N_KERNELS_HPP

#include<cassert>
#include <complex>

namespace kernels 
{

void fill_n(int * first, int N, int stride, int const value);
void fill_n(float * first, int N,  int stride, float const value);
void fill_n(double * first, int N,  int stride, double const value);
void fill_n(std::complex<float> * first, int N,  int stride, std::complex<float> const value);
void fill_n(std::complex<double> * first, int N,  int stride, std::complex<double> const value);

void fill_n(int * first, int N, int const value);
void fill_n(float * first, int N, float const value);
void fill_n(double * first, int N, double const value);
void fill_n(std::complex<float> * first, int N, std::complex<float> const value);
void fill_n(std::complex<double> * first, int N, std::complex<double> const value);


}

#endif
