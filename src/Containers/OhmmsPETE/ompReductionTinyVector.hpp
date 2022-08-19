//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_OMPREDUCTION_TINYVECTOR_H
#define QMCPLUSPLUS_OMPREDUCTION_TINYVECTOR_H

#include <complex>
#include "config.h"
#include "TinyVector.h"

namespace qmcplusplus
{
#if !defined(OPENMP_NO_UDR)
PRAGMA_OFFLOAD("omp declare reduction(+: TinyVector<float, OHMMS_DIM>: omp_out += omp_in)")
PRAGMA_OFFLOAD("omp declare reduction(+: TinyVector<double, OHMMS_DIM>: omp_out += omp_in)")
#endif

#if !defined(OPENMP_NO_COMPLEX) && !defined(OPENMP_NO_UDR)
PRAGMA_OFFLOAD("omp declare reduction(+: TinyVector<std::complex<float>, OHMMS_DIM>: omp_out += omp_in)")
PRAGMA_OFFLOAD("omp declare reduction(+: TinyVector<std::complex<double>, OHMMS_DIM>: omp_out += omp_in)")
#endif
}
#endif // QMCPLUSPLUS_OMPREDUCTION_TINYVECTOR_H
