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

#ifndef AFQMC_GET_DIAGONAL_KERNELS_HPP
#define AFQMC_GET_DIAGONAL_KERNELS_HPP

#include <cassert>
#include <complex>

namespace kernels
{
void get_diagonal_strided(int nk,
                          int ni,
                          std::complex<double> const* B,
                          int ldb,
                          int stride,
                          std::complex<double>* A,
                          int lda);
void get_diagonal_strided(int nk,
                          int ni,
                          std::complex<float> const* B,
                          int ldb,
                          int stride,
                          std::complex<float>* A,
                          int lda);


} // namespace kernels

#endif
