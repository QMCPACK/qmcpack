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

#ifndef AFQMC_AJW_TO_WAJ_H
#define AFQMC_AJW_TO_WAJ_H

#include <complex>

namespace kernels
{
void ajw_to_waj(int na, int nj, int nw, int inca, double const* A, double* B);
void ajw_to_waj(int na, int nj, int nw, int inca, float const* A, float* B);
void ajw_to_waj(int na, int nj, int nw, int inca, std::complex<double> const* A, std::complex<double>* B);
void ajw_to_waj(int na, int nj, int nw, int inca, std::complex<float> const* A, std::complex<float>* B);


void transpose_wabn_to_wban(int nwalk,
                            int na,
                            int nb,
                            int nchol,
                            std::complex<double> const* Tab,
                            std::complex<double>* Tba);
void transpose_wabn_to_wban(int nwalk,
                            int na,
                            int nb,
                            int nchol,
                            std::complex<float> const* Tab,
                            std::complex<float>* Tba);
void transpose_wabn_to_wban(int nwalk,
                            int na,
                            int nb,
                            int nchol,
                            std::complex<double> const* Tab,
                            std::complex<float>* Tba);
void transpose_wabn_to_wban(int nwalk,
                            int na,
                            int nb,
                            int nchol,
                            std::complex<float> const* Tab,
                            std::complex<double>* Tba);

} // namespace kernels
#endif
