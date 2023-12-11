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

#ifndef AFQMC_DOT_WABN_H
#define AFQMC_DOT_WABN_H

#include <cassert>
#include <complex>
#include "AFQMC/Numerics/detail/HIP/Kernels/hip_settings.h"

namespace kernels
{
void dot_wabn(int nwalk,
              int nocc,
              int nchol,
              std::complex<double> const alpha,
              std::complex<double> const* Tab,
              std::complex<double>* y,
              int incy);
void dot_wabn(int nwalk,
              int nocc,
              int nchol,
              std::complex<float> const alpha,
              std::complex<float> const* Tab,
              std::complex<float>* y,
              int incy);
void dot_wabn(int nwalk,
              int nocc,
              int nchol,
              std::complex<float> const alpha,
              std::complex<float> const* Tab,
              std::complex<double>* y,
              int incy);


void dot_wanb(int nwalk,
              int nocc,
              int nchol,
              std::complex<double> const alpha,
              std::complex<double> const* Tab,
              std::complex<double>* y,
              int incy);
void dot_wanb(int nwalk,
              int nocc,
              int nchol,
              std::complex<float> const alpha,
              std::complex<float> const* Tab,
              std::complex<float>* y,
              int incy);
void dot_wanb(int nwalk,
              int nocc,
              int nchol,
              std::complex<float> const alpha,
              std::complex<float> const* Tab,
              std::complex<double>* y,
              int incy);


void dot_wpan_waqn_Fwpq(int nwalk,
                        int nmo,
                        int nchol,
                        std::complex<double> const alpha,
                        std::complex<double> const* Tab,
                        std::complex<double>* F);
void dot_wpan_waqn_Fwpq(int nwalk,
                        int nmo,
                        int nchol,
                        std::complex<float> const alpha,
                        std::complex<float> const* Tab,
                        std::complex<double>* F);
void dot_wpan_waqn_Fwpq(int nwalk,
                        int nmo,
                        int nchol,
                        std::complex<float> const alpha,
                        std::complex<float> const* Tab,
                        std::complex<float>* F);

} // namespace kernels

#endif
