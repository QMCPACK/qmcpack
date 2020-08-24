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

#ifndef AFQMC_ADIAGAPY_KERNELS_HPP
#define AFQMC_ADIAGAPY_KERNELS_HPP

#include <cassert>
#include <complex>

namespace kernels
{
void adiagApy(int N, float const alpha, float const* A, int lda, float* y, int incy);
void adiagApy(int N,
              std::complex<float> const alpha,
              std::complex<float> const* A,
              int lda,
              std::complex<float>* y,
              int incy);
void adiagApy(int N, double const alpha, double const* A, int lda, double* y, int incy);
void adiagApy(int N,
              std::complex<double> const alpha,
              std::complex<double> const* A,
              int lda,
              std::complex<double>* y,
              int incy);

} // namespace kernels

#endif
