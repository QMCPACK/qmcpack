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

#ifndef AFQMC_FILL_N_KERNELS_HPP
#define AFQMC_FILL_N_KERNELS_HPP

#include <cassert>
#include <complex>

namespace kernels
{
void fill_n(char* first, int N, int stride, char const value);
void fill_n(int* first, int N, int stride, int const value);
void fill_n(float* first, int N, int stride, float const value);
void fill_n(double* first, int N, int stride, double const value);
void fill_n(std::complex<float>* first, int N, int stride, std::complex<float> const value);
void fill_n(std::complex<double>* first, int N, int stride, std::complex<double> const value);

void fill_n(char* first, int N, char const value);
void fill_n(int* first, int N, int const value);
void fill_n(long unsigned int* first, long unsigned int N, const long unsigned int value);
void fill_n(float* first, int N, float const value);
void fill_n(double* first, int N, double const value);
void fill_n(std::complex<float>* first, int N, std::complex<float> const value);
void fill_n(std::complex<double>* first, int N, std::complex<double> const value);

void fill2D_n(int N, int M, int* A, int lda, int const value);
void fill2D_n(int N, int M, float* A, int lda, float const value);
void fill2D_n(int N, int M, double* A, int lda, double const value);
void fill2D_n(int N, int M, std::complex<double>* A, int lda, std::complex<double> const value);
void fill2D_n(int N, int M, std::complex<float>* A, int lda, std::complex<float> const value);

} // namespace kernels

#endif
