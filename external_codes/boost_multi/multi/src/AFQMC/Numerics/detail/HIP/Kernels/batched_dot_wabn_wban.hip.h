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

#ifndef AFQMC_BATCHED_DOT_WABN_WBAN_H
#define AFQMC_BATCHED_DOT_WABN_WBAN_H

#include <cassert>
#include <complex>
#include "AFQMC/Numerics/detail/HIP/Kernels/hip_settings.h"

namespace kernels
{
void batched_dot_wabn_wban(int nbatch,
                           int nwalk,
                           int nocc,
                           int nchol,
                           std::complex<double> const* alpha,
                           std::complex<double> const* Tab,
                           std::complex<double>* y,
                           int incy);
void batched_dot_wabn_wban(int nbatch,
                           int nwalk,
                           int nocc,
                           int nchol,
                           std::complex<float> const* alpha,
                           std::complex<float> const* Tab,
                           std::complex<float>* y,
                           int incy);
void batched_dot_wabn_wban(int nbatch,
                           int nwalk,
                           int nocc,
                           int nchol,
                           std::complex<float> const* alpha,
                           std::complex<float> const* Tab,
                           std::complex<double>* y,
                           int incy);

void batched_dot_wanb_wbna(int nbatch,
                           int nwalk,
                           int nocc,
                           int nchol,
                           std::complex<double> const* alpha,
                           std::complex<double> const* Tab,
                           std::complex<double>* y,
                           int incy);
void batched_dot_wanb_wbna(int nbatch,
                           int nwalk,
                           int nocc,
                           int nchol,
                           std::complex<float> const* alpha,
                           std::complex<float> const* Tab,
                           std::complex<float>* y,
                           int incy);
void batched_dot_wanb_wbna(int nbatch,
                           int nwalk,
                           int nocc,
                           int nchol,
                           std::complex<float> const* alpha,
                           std::complex<float> const* Tab,
                           std::complex<double>* y,
                           int incy);

} // namespace kernels

#endif
