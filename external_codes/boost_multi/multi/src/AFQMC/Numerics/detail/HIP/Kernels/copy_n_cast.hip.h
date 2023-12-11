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

#ifndef AFQMC_COPY_N_CAST_KERNELS_HPP
#define AFQMC_COPY_N_CAST_KERNELS_HPP

#include <stdexcept>
#include <cassert>
#include <complex>
#include <iostream>

namespace kernels
{
void copy_n_cast(double const* A, int n, float* B);
void copy_n_cast(float const* A, int n, double* B);
void copy_n_cast(std::complex<double> const* A, int n, std::complex<float>* B);
void copy_n_cast(std::complex<float> const* A, int n, std::complex<double>* B);
inline void copy_n_cast(std::complex<float> const* A, int n, std::complex<float>* B)
{
  std::cerr << " Should not be calling copy_n_cast<T,T>. \n" << std::endl;
  throw std::runtime_error("Calling cast_n_copy(float const*,n,float*)");
}
inline void copy_n_cast(std::complex<double> const* A, int n, std::complex<double>* B)
{
  std::cerr << " Should not be calling copy_n_cast<T,T>. \n" << std::endl;
  throw std::runtime_error("Calling cast_n_copy(double const*,n,double*)");
}

} // namespace kernels

#endif
