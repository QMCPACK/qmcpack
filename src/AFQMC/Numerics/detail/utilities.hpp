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

#include <complex>
#include "AFQMC/config.0.h"
#include "AFQMC/Memory/raw_pointers.hpp"
#include "AFQMC/Memory/SharedMemory/shm_ptr_with_raw_ptr_dispatch.hpp"

namespace ma
{
using qmcplusplus::afqmc::to_address;

enum TENSOR_OPERATIONS
{
  TOp_PLUS,
  TOp_MINUS,
  TOp_MUL,
  TOp_DIV
};

static const int INCX                   = 1;
static const int INCY                   = 1;
static const char UPLO                  = 'L';
static const char TRANS                 = 'T';
static const char NOTRANS               = 'N';
static const float sone                 = 1.0;
static const float szero                = 0.0;
static const double done                = 1.0;
static const double dzero               = 0.0;
static const std::complex<float> cone   = std::complex<float>(1.0, 0.0);
static const std::complex<float> czero  = std::complex<float>(0.0, 0.0);
static const std::complex<double> zone  = std::complex<double>(1.0, 0.0);
static const std::complex<double> zzero = std::complex<double>(0.0, 0.0);

/*
#if defined(HAVE_MKL)
using CBLAS_LAYOUT = enum {CblasRowMajor=101, CblasColMajor=102};
using CBLAS_TRANSPOSE = enum {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
#endif
*/

inline double const& real(double const& d) { return d; }
inline float const& real(float const& f) { return f; }

inline double conj(double const& d) { return d; }
inline float conj(float const& f) { return f; }

inline std::complex<double> conj(std::complex<double> const& d) { return std::conj(d); }
inline std::complex<float> conj(std::complex<float> const& f) { return std::conj(f); }
//template<typename T>
//T conj(T const& v) { return v; }
//template<typename T>
//std::complex<T> conj(std::complex<T> const& v) { return std::conj(v); }

template<class Ptr>
auto pointer_dispatch(Ptr p)
{
  return p;
}

template<typename T>
T* pointer_dispatch(shm::shm_ptr_with_raw_ptr_dispatch<T> p)
{
  return to_address(p);
}


} // namespace ma

#endif
