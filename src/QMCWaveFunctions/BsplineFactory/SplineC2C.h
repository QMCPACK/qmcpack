//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////


/** @file
 *
 * class to handle complex splines to complex orbitals with splines of arbitrary precision
 */
#ifndef QMCPLUSPLUS_SPLINE_C2C_H
#define QMCPLUSPLUS_SPLINE_C2C_H

#include "Configuration.h"
#include "QMCWaveFunctions/BsplineFactory/SplineC2CT.h"

namespace qmcplusplus
{
template<class ST>
using SplineC2C = SplineC2CT<ST, QMCTraits::ValueType>;

} // namespace qmcplusplus
#endif
