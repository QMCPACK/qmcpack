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


/** @file DynamicCoordinatesBuilder.h
 */
#ifndef QMCPLUSPLUS_QUANTUM_VARIABLE_BUILDER_H
#define QMCPLUSPLUS_QUANTUM_VARIABLE_BUILDER_H

#include "Particle/DynamicCoordinates.h"

namespace qmcplusplus
{
/** create DynamicCoordinates based on kind
 */
std::unique_ptr<DynamicCoordinates> createDynamicCoordinates(
    const DynamicCoordinateKind kind = DynamicCoordinateKind::DC_POS);
} // namespace qmcplusplus
#endif
