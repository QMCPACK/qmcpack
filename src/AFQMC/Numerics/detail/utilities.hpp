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

#ifndef AFQMC_MA_UTILITIES_HPP
#define AFQMC_MA_UTILITIES_HPP

#include<complex>
#include "AFQMC/Memory/raw_pointers.hpp"

namespace ma {

using qmcplusplus::afqmc::to_address;

inline double const& real(double const& d){return d;}
inline float const& real(float const& f){return f;}

inline double conj(double const& d){return d;}
inline float conj(float const& f){return f;}

inline std::complex<double> conj(std::complex<double> const& d){return std::conj(d);}
inline std::complex<float>  conj(std::complex<float> const& f){return std::conj(f);}
}

#endif
