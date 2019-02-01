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

#ifndef AFQMC_BATCHEDDOT_KERNELS_HPP
#define AFQMC_BATCHEDDOT_KERNELS_HPP

#include <complex>

namespace kernels
{

void ajw_to_waj(int na, int nj, int nw, int inca,  
                double const* A, 
                double * B); 
void ajw_to_waj(int na, int nj, int nw, int inca,
                float const* A,
                float * B);
void ajw_to_waj(int na, int nj, int nw, int inca,
                std::complex<double> const* A,
                std::complex<double> * B);
void ajw_to_waj(int na, int nj, int nw, int inca,
                std::complex<float> const* A,
                std::complex<float> * B);

}
#endif
