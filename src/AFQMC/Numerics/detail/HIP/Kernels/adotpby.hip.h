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

#ifndef AFQMC_ADOTPBY_KERNELS_HPP
#define AFQMC_ADOTPBY_KERNELS_HPP

#include <cassert>
#include <complex>

namespace kernels
{
void adotpby(int N,
             double const alpha,
             double const* x,
             int const incx,
             double const* y,
             int const incy,
             double const beta,
             double* res);

void adotpby(int N,
             std::complex<double> const alpha,
             std::complex<double> const* x,
             int const incx,
             std::complex<double> const* y,
             int const incy,
             std::complex<double> const beta,
             std::complex<double>* res);

void adotpby(int N,
             float const alpha,
             float const* x,
             int const incx,
             float const* y,
             int const incy,
             float const beta,
             float* res);

void adotpby(int N,
             std::complex<float> const alpha,
             std::complex<float> const* x,
             int const incx,
             std::complex<float> const* y,
             int const incy,
             std::complex<float> const beta,
             std::complex<float>* res);

void adotpby(int N,
             float const alpha,
             float const* x,
             int const incx,
             float const* y,
             int const incy,
             double const beta,
             double* res);

void adotpby(int N,
             std::complex<float> const alpha,
             std::complex<float> const* x,
             int const incx,
             std::complex<float> const* y,
             int const incy,
             std::complex<double> const beta,
             std::complex<double>* res);

void strided_adotpby(int NB,
                     int N,
                     std::complex<double> const alpha,
                     std::complex<double> const* A,
                     int const lda,
                     std::complex<double> const* B,
                     int const ldb,
                     std::complex<double> const beta,
                     std::complex<double>* C,
                     int ldc);

void strided_adotpby(int NB,
                     int N,
                     std::complex<float> const alpha,
                     std::complex<float> const* A,
                     int const lda,
                     std::complex<float> const* B,
                     int const ldb,
                     std::complex<double> const beta,
                     std::complex<double>* C,
                     int ldc);


} // namespace kernels

#endif
