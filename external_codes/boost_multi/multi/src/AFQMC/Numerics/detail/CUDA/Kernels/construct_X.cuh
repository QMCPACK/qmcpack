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

#ifndef AFQMC_CONSTRUCT_X_H
#define AFQMC_CONSTRUCT_X_H

#include <cassert>
#include <complex>

namespace kernels
{
void construct_X(int nCV,
                 int nsteps,
                 int nwalk,
                 bool free_projection,
                 double sqrtdt,
                 double vbound,
                 std::complex<double> const* vMF,
                 std::complex<double> const* vbias,
                 std::complex<double>* HW,
                 std::complex<double>* MF,
                 std::complex<double>* X);
void construct_X(int nCV,
                 int nsteps,
                 int nwalk,
                 bool free_projection,
                 double sqrtdt,
                 double vbound,
                 std::complex<double> const* vMF,
                 std::complex<float> const* vbias,
                 std::complex<double>* HW,
                 std::complex<double>* MF,
                 std::complex<float>* X);
void construct_X(int nCV,
                 int nsteps,
                 int nwalk,
                 bool free_projection,
                 double sqrtdt,
                 double vbound,
                 std::complex<double> const* vMF,
                 std::complex<float> const* vbias,
                 std::complex<double>* HW,
                 std::complex<double>* MF,
                 std::complex<double>* X);
void construct_X(int nCV,
                 int nsteps,
                 int nwalk,
                 bool free_projection,
                 double sqrtdt,
                 double vbound,
                 std::complex<double> const* vMF,
                 std::complex<double> const* vbias,
                 std::complex<double>* HW,
                 std::complex<double>* MF,
                 std::complex<float>* X);


} // namespace kernels

#endif
