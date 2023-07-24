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


/** @file RealSpacePostions.h
 */
#ifndef QMCPLUSPLUS_REALSPACE_POSITIONS_H
#define QMCPLUSPLUS_REALSPACE_POSITIONS_H

#include "Configuration.h"
#include "Particle/RealSpacePositionsT.h"

namespace qmcplusplus
{
using RealSpacePositions = RealSpacePositionsT<QMCTraits::ValueType>;

} // namespace qmcplusplus
#endif
