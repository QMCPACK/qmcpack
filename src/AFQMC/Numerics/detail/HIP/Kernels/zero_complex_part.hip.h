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

#ifndef AFQMC_ZERO_COMPLEX_PART_KERNELS_HPP
#define AFQMC_ZERO_COMPLEX_PART_KERNELS_HPP

#include <complex>

namespace kernels
{
void zero_complex_part(int n, std::complex<double>* x);
void zero_complex_part(int n, std::complex<float>* x);
void zero_complex_part(int n, double* x);
void zero_complex_part(int n, float* x);

} // namespace kernels

#endif
