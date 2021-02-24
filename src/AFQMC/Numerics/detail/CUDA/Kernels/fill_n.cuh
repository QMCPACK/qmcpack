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

#ifndef AFQMC_FILL_N_KERNELS_HPP
#define AFQMC_FILL_N_KERNELS_HPP

#include <cassert>
#include <complex>

namespace kernels
{

template<typename T>
void fill_n(T* first, long N, long stride, T const value);

template<typename T>
void fill_n(T* first, long N, T const value);

template<typename T>
void fill2D_n(long N, long M, T* A, long lda, T const value);

} // namespace kernels

#endif
