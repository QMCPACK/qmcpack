//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SIMULATIONCELL_H
#define QMCPLUSPLUS_SIMULATIONCELL_H

#include "Configuration.h"
#include "SimulationCellT.h"

namespace qmcplusplus
{
using SimulationCell = SimulationCellT<QMCTraits::ValueType>;

} // namespace qmcplusplus
#endif
