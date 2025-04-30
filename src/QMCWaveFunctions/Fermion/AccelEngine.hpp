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

#ifndef QMCPLUSPLUS_ACCELENGINE_HPP
#define QMCPLUSPLUS_ACCELENGINE_HPP

#include "config.h"
#include "QMCWaveFunctions/Fermion/DelayedUpdate.h"
#if defined(ENABLE_CUDA) || defined(ENABLE_SYCL)
#include "QMCWaveFunctions/Fermion/DelayedUpdateAccel.h"
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
#endif

namespace qmcplusplus
{
template<PlatformKind P, typename T, typename FP_T>
struct AccelEngine;

template<typename T, typename FP_T>
struct AccelEngine<PlatformKind::CPU, T, FP_T>
{
  static constexpr bool inverter_supported = false;
  DelayedUpdate<T> update_eng_;
};

#if defined(ENABLE_CUDA)
template<typename T, typename FP_T>
struct AccelEngine<PlatformKind::CUDA, T, FP_T>
{
  static constexpr bool inverter_supported = true;
  DelayedUpdateAccel<PlatformKind::CUDA, T> update_eng_;
#if defined(QMC_CUDA2HIP)
  using Inverter = rocSolverInverter<FP_T>;
#else
  using Inverter = cuSolverInverter<FP_T>;
#endif
  Inverter inverter_;
};
#endif

#if defined(ENABLE_SYCL)
template<typename T, typename FP_T>
struct AccelEngine<PlatformKind::SYCL, T, FP_T>
{
  static constexpr bool inverter_supported = true;
  DelayedUpdateAccel<PlatformKind::SYCL, T> update_eng_;
  using Inverter = syclSolverInverter<FP_T>;
  Inverter inverter_;
};
#endif

} // namespace qmcplusplus

#endif
