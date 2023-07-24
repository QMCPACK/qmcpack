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


/** @file HybridRepReal.h
 *
 * hold HybridRepReal
 */
#ifndef QMCPLUSPLUS_HYBRIDREP_REAL_H
#define QMCPLUSPLUS_HYBRIDREP_REAL_H

#include "Configuration.h"
#include "QMCWaveFunctions/BsplineFactory/HybridRepRealT.h"

namespace qmcplusplus
{
template<typename SPLINEBASE>
using HybridRepReal = HybridRepRealT<SPLINEBASE>;

} // namespace qmcplusplus
#endif
