//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Youngjun Lee, leey@anl.gov, Argonne National Laboratory
//
// File created by: Youngjun Lee, leey@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_INVERTERACCEL_HPP
#define QMCPLUSPLUS_INVERTERACCEL_HPP

#include "config.h"
#if defined(ENABLE_CUDA)
#if defined(QMC_CUDA2HIP)
#include "rocSolverInverter.hpp"
#else
#include "cuSolverInverter.hpp"
#endif
#endif
#if defined(ENABLE_SYCL)
#include "syclSolverInverter.hpp"
#endif

namespace qmcplusplus
{
template<PlatformKind P, typename T>
struct InverterAccel;

#if defined(ENABLE_CUDA)
template<typename T>
struct InverterAccel<PlatformKind::CUDA, T>
{
#if defined(QMC_CUDA2HIP)
  using Inverter = rocSolverInverter<T>;
#else
  using Inverter = cuSolverInverter<T>;
#endif
};
#endif

#if defined(ENABLE_SYCL)
template<typename T>
struct InverterAccel<PlatformKind::SYCL, T>
{
  using Inverter = syclSolverInverter<T>;
};
#endif

} // namespace qmcplusplus

#endif
