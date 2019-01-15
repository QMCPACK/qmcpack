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

#ifndef AFQMC_DOT_WABN_WBAN_H
#define AFQMC_DOT_WABN_WBAN_H

#include<cassert>
#include <complex>

namespace kernels 
{

void dot_wabn_wban( int nw, int na, int nb, int nc,
                    double const alpha, double const* A, double const* B, double* y, int incy);
void dot_wabn_wban( int nw, int na, int nb, int nc,
                    double const alpha, float const* A, float const* B, double* y, int incy);
void dot_wabn_wban( int nw, int na, int nb, int nc,
                    std::complex<double> const alpha, std::complex<double> const* A, 
                    std::complex<double> const* B, std::complex<double>* y, int incy);
void dot_wabn_wban( int nw, int na, int nb, int nc,
                    std::complex<double> const alpha, std::complex<float> const* A, 
                    std::complex<float> const* B, std::complex<double>* y, int incy);
}

#endif
