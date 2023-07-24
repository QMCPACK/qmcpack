//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file RealSpacePostionsOMPTarget.h
 */
#ifndef QMCPLUSPLUS_REALSPACE_POSITIONS_OMPTARGET_H
#define QMCPLUSPLUS_REALSPACE_POSITIONS_OMPTARGET_H

#include "Configuration.h"
#include "Particle/RealSpacePositionsTOMPTarget.h"

namespace qmcplusplus
{
using RealSpacePositionsOMPTarget = RealSpacePositionsTOMPTarget<QMCTraits::ValueType>;

} // namespace qmcplusplus
#endif
