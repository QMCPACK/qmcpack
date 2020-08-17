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

#ifndef AFQMC_ACAXPBB_KERNELS_HPP
#define AFQMC_ACAXPBB_KERNELS_HPP

#include <cassert>
#include <complex>

namespace kernels
{
void acAxpbB(int m,
             int n,
             double const alpha,
             double const* A,
             int lda,
             double const* x,
             int incx,
             double const beta,
             double* B,
             int ldb);

void acAxpbB(int m,
             int n,
             float const alpha,
             float const* A,
             int lda,
             float const* x,
             int incx,
             float const beta,
             float* B,
             int ldb);

void acAxpbB(int m,
             int n,
             std::complex<double> const alpha,
             std::complex<double> const* A,
             int lda,
             std::complex<double> const* x,
             int incx,
             std::complex<double> const beta,
             std::complex<double>* B,
             int ldb);

void acAxpbB(int m,
             int n,
             std::complex<float> const alpha,
             std::complex<float> const* A,
             int lda,
             std::complex<float> const* x,
             int incx,
             std::complex<float> const beta,
             std::complex<float>* B,
             int ldb);

} // namespace kernels

#endif
