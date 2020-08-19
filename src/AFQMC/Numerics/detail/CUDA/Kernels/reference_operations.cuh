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

#include <complex>

namespace kernels
{
// +=
void op_plus(double* x, double inc);
void op_plus(float* x, float inc);
void op_plus(std::complex<double>* x, std::complex<double> inc);
void op_plus(std::complex<float>* x, std::complex<float> inc);

// -=
void op_minus(double* x, double inc);
void op_minus(float* x, float inc);
void op_minus(std::complex<double>* x, std::complex<double> inc);
void op_minus(std::complex<float>* x, std::complex<float> inc);

// *=
void op_times(double* x, double inc);
void op_times(float* x, float inc);
void op_times(std::complex<double>* x, std::complex<double> inc);
void op_times(std::complex<float>* x, std::complex<float> inc);

// /=
void op_div(double* x, double inc);
void op_div(float* x, float inc);
void op_div(std::complex<double>* x, std::complex<double> inc);
void op_div(std::complex<float>* x, std::complex<float> inc);

} // namespace kernels
