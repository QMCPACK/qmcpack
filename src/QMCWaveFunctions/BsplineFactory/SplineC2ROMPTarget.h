//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file SplineC2ROMPTarget.h
 *
 * class to handle complex splines to real orbitals with splines of arbitrary precision
 * splines storage and computation is offloaded to accelerators using OpenMP target
 */
#ifndef QMCPLUSPLUS_SPLINE_C2R_OMPTARGET_H
#define QMCPLUSPLUS_SPLINE_C2R_OMPTARGET_H

#include "QMCWaveFunctions/BsplineFactory/SplineC2ROMPTargetT.h"

namespace qmcplusplus
{
template<class ST>
using SplineC2ROMPTarget = SplineC2ROMPTargetT<ST, QMCTraits::ValueType>;

} // namespace qmcplusplus
#endif
