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

#ifndef AFQMC_INPLACE_PRODUCT_KERNELS_HPP
#define AFQMC_INPLACE_PRODUCT_KERNELS_HPP

#include <cassert>
#include <complex>

namespace kernels
{
void inplace_product(int nbatch, int n, int m, double const* B, int ldb, std::complex<double>* A, int lda);
void inplace_product(int nbatch, int n, int m, float const* B, int ldb, std::complex<float>* A, int lda);
void inplace_product(int nbatch,
                     int n,
                     int m,
                     std::complex<double> const* B,
                     int ldb,
                     std::complex<double>* A,
                     int lda);
void inplace_product(int nbatch, int n, int m, std::complex<float> const* B, int ldb, std::complex<float>* A, int lda);

} // namespace kernels

#endif
