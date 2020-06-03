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

#ifndef AFQMC_TAB_TO_KL_H
#define AFQMC_TAB_TO_KL_H

#include<cassert>
#include <complex>

namespace kernels 
{

void Tab_to_Kl(int nwalk, int nocc, int nchol, std::complex<float> const* Tab,
                    std::complex<float> * Kl);

void Tab_to_Kl(int nwalk, int nocc, int nchol, std::complex<double> const* Tab,
                    std::complex<double> * Kl);

void Tanb_to_Kl(int nwalk, int nocc, int nchol, int nchol_tot, std::complex<float> const* Tab,
                    std::complex<float> * Kl);

void Tanb_to_Kl(int nwalk, int nocc, int nchol, int nchol_tot, std::complex<double> const* Tab,
                    std::complex<double> * Kl);

}

#endif
