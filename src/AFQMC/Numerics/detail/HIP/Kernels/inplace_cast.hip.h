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

#ifndef AFQMC_INPLACE_CAST_KERNELS_HPP
#define AFQMC_INPLACE_CAST_KERNELS_HPP

#include <cassert>
#include <complex>

namespace kernels
{
void inplace_cast(unsigned long n, std::complex<float>* A, std::complex<double>* B);
void inplace_cast(unsigned long n, std::complex<double>* A, std::complex<float>* B);
void inplace_cast(long n, std::complex<float>* A, std::complex<double>* B);
void inplace_cast(long n, std::complex<double>* A, std::complex<float>* B);

} // namespace kernels

#endif
