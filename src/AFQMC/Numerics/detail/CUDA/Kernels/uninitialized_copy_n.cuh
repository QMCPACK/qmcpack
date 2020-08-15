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

#ifndef AFQMC_UNINITIALIZED_COPY_N_KERNELS_HPP
#define AFQMC_UNINITIALIZED_COPY_N_KERNELS_HPP

#include <cassert>
#include <complex>

namespace kernels
{
void uninitialized_copy_n(int N, double const* first, int incx, double* array, int incy);
void uninitialized_copy_n(int N, std::complex<double> const* first, int incx, std::complex<double>* array, int incy);
void uninitialized_copy_n(int N, int const* first, int incx, int* array, int incy);

// long
void uninitialized_copy_n(long N, double const* first, long incx, double* array, long incy);
void uninitialized_copy_n(long N, std::complex<double> const* first, long incx, std::complex<double>* array, long incy);
void uninitialized_copy_n(long N, int const* first, long incx, int* array, long incy);

} // namespace kernels

#endif
