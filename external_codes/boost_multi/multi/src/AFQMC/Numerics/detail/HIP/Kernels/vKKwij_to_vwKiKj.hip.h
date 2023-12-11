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

#ifndef AFQMC_VKKWIJ_TO_VWKIKJ_H
#define AFQMC_VKKWIJ_TO_VWKIKJ_H

#include <complex>

namespace kernels
{
void vKKwij_to_vwKiKj(int nwalk,
                      int nkpts,
                      int nmo_max,
                      int nmo_tot,
                      int* kk,
                      int* nmo,
                      int* nmo0,
                      double const* A,
                      double* B);
void vKKwij_to_vwKiKj(int nwalk,
                      int nkpts,
                      int nmo_max,
                      int nmo_tot,
                      int* kk,
                      int* nmo,
                      int* nmo0,
                      float const* A,
                      float* B);
void vKKwij_to_vwKiKj(int nwalk,
                      int nkpts,
                      int nmo_max,
                      int nmo_tot,
                      int* kk,
                      int* nmo,
                      int* nmo0,
                      float const* A,
                      double* B);
void vKKwij_to_vwKiKj(int nwalk,
                      int nkpts,
                      int nmo_max,
                      int nmo_tot,
                      int* kk,
                      int* nmo,
                      int* nmo0,
                      std::complex<double> const* A,
                      std::complex<double>* B);
void vKKwij_to_vwKiKj(int nwalk,
                      int nkpts,
                      int nmo_max,
                      int nmo_tot,
                      int* kk,
                      int* nmo,
                      int* nmo0,
                      std::complex<float> const* A,
                      std::complex<float>* B);
void vKKwij_to_vwKiKj(int nwalk,
                      int nkpts,
                      int nmo_max,
                      int nmo_tot,
                      int* kk,
                      int* nmo,
                      int* nmo0,
                      std::complex<float> const* A,
                      std::complex<double>* B);
} // namespace kernels

#endif
