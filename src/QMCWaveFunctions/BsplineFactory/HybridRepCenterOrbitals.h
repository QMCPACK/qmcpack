//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////

/** @file HybridRepCenterOrbitals.h
 *
 * Hybrid representation atomic centered orbitals
 */
#ifndef QMCPLUSPLUS_HYBRIDREP_CENTER_ORBITALS_H
#define QMCPLUSPLUS_HYBRIDREP_CENTER_ORBITALS_H

#include "Configuration.h"
#include "QMCWaveFunctions/BsplineFactory/HybridRepCenterOrbitalsT.h"

namespace qmcplusplus
{
template<typename ST>
using AtomicOrbitals = AtomicOrbitalsT<ST>;

template<typename ST>
using HybridRepCenterOrbitals = HybridRepCenterOrbitals<ST>;
} // namespace qmcplusplus
#endif
