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

#ifndef AFQMC_BATCHED_TAB_TO_KLR_H
#define AFQMC_BATCHED_TAB_TO_KLR_H

#include <cassert>
#include <complex>

namespace kernels
{
void batched_Tab_to_Klr(int nterms,
                        int nwalk,
                        int nocc,
                        int nchol_max,
                        int nchol_tot,
                        int ncholQ,
                        int ncholQ0,
                        int* kdiag,
                        std::complex<float> const* Tab,
                        std::complex<float>* Kl,
                        std::complex<float>* Kr);

void batched_Tab_to_Klr(int nterms,
                        int nwalk,
                        int nocc,
                        int nchol_max,
                        int nchol_tot,
                        int ncholQ,
                        int ncholQ0,
                        int* kdiag,
                        std::complex<double> const* Tab,
                        std::complex<double>* Kl,
                        std::complex<double>* Kr);

void batched_Tanb_to_Klr(int nterms,
                         int nwalk,
                         int nocc,
                         int nchol_max,
                         int nchol_tot,
                         int ncholQ,
                         int ncholQ0,
                         int* kdiag,
                         std::complex<float> const* Tab,
                         std::complex<float>* Kl,
                         std::complex<float>* Kr);

void batched_Tanb_to_Klr(int nterms,
                         int nwalk,
                         int nocc,
                         int nchol_max,
                         int nchol_tot,
                         int ncholQ,
                         int ncholQ0,
                         int* kdiag,
                         std::complex<double> const* Tab,
                         std::complex<double>* Kl,
                         std::complex<double>* Kr);

} // namespace kernels

#endif
