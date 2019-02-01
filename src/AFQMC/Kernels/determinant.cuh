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

#ifndef AFQMC_DETERMINANT_KERNELS_HPP
#define AFQMC_DETERMINANT_KERNELS_HPP

#include<cassert>
#include <complex>

namespace kernels 
{

//double determinant_from_getrf_gpu(int N, double *m, int lda, int *piv);
//std::complex<double> determinant_from_getrf_gpu(int N, std::complex<double> *m, int lda, int *piv);

void determinant_from_getrf_gpu(int N, double *m, int lda, int *piv, double* res);
void determinant_from_getrf_gpu(int N, std::complex<double> *m, int lda, int *piv, std::complex<double>* res);


}

#endif
