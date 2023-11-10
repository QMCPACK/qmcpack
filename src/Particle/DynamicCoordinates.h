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


/** @file DynamicCoordinates.h
 */
#ifndef QMCPLUSPLUS_DYNAMICCOORDINATES_H
#define QMCPLUSPLUS_DYNAMICCOORDINATES_H

#include "Configuration.h"
#include "Particle/DynamicCoordinatesT.h"

namespace qmcplusplus
{
using DynamicCoordinates = DynamicCoordinatesT<QMCTraits::ValueType>;
} // namespace qmcplusplus
#endif
