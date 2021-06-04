//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_OMPREDUCTION_H
#define QMCPLUSPLUS_OMPREDUCTION_H

#include <complex>
#include "config.h"

#if !defined(OPENMP_NO_COMPLEX) && !defined(OPENMP_NO_UDR)
PRAGMA_OFFLOAD("omp declare reduction(+: std::complex<float>: omp_out += omp_in)")
PRAGMA_OFFLOAD("omp declare reduction(+: std::complex<double>: omp_out += omp_in)")
#endif

#endif // QMCPLUSPLUS_OMPREDUCTION_H
