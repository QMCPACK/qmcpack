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

#ifndef AFQMC_PRINT_KERNELS_HPP
#define AFQMC_PRINT_KERNELS_HPP

#include <cassert>
#include <complex>

namespace kernels
{
void print(std::string str, std::complex<double> const* p, int n);
void print(std::string str, double const* p, int n);
void print(std::string str, int const* p, int n);
void print(std::string str, size_t const* p, int n);
void print(std::string str, long const* p, int n);


} // namespace kernels

#endif
