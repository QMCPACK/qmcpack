//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@intel.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Anouar Benali, benali@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file
 *
 * class to handle complex splines to real orbitals with splines of arbitrary precision
 */
#ifndef QMCPLUSPLUS_SPLINE_C2R_H
#define QMCPLUSPLUS_SPLINE_C2R_H

#include "QMCWaveFunctions/BsplineFactory/SplineC2RT.h"

namespace qmcplusplus
{
template<class ST>
using SplineC2R = SplineC2RT<ST, QMCTraits::ValueType>;

} // namespace qmcplusplus
#endif
