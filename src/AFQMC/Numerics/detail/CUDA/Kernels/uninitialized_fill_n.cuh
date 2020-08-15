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

#ifndef AFQMC_UNINITIALIZED_FILL_N_KERNELS_HPP
#define AFQMC_UNINITIALIZED_FILL_N_KERNELS_HPP

#include <cassert>
#include <complex>

namespace kernels
{
void uninitialized_fill_n(bool* first, int N, bool const value);
void uninitialized_fill_n(int* first, int N, int const value);
void uninitialized_fill_n(float* first, int N, float const value);
void uninitialized_fill_n(double* first, int N, double const value);
void uninitialized_fill_n(std::complex<float>* first, int N, std::complex<float> const value);
void uninitialized_fill_n(std::complex<double>* first, int N, std::complex<double> const value);
//void uninitialized_fill_n(double2 * first, int N,  double2 const value);

void uninitialized_fill_n(bool* first, long N, bool const value);
void uninitialized_fill_n(int* first, long N, int const value);
void uninitialized_fill_n(float* first, long N, float const value);
void uninitialized_fill_n(double* first, long N, double const value);
void uninitialized_fill_n(std::complex<float>* first, long N, std::complex<float> const value);
void uninitialized_fill_n(std::complex<double>* first, long N, std::complex<double> const value);
//void uninitialized_fill_n(double2 * first, long N,  double2 const value);

} // namespace kernels

#endif
