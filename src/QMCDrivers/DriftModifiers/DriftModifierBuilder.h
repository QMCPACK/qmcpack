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


#ifndef QMCPLUSPLUS_DRIFTMODIFIER_BUILDER_H
#define QMCPLUSPLUS_DRIFTMODIFIER_BUILDER_H

#include "QMCDrivers/DriftModifiers/DriftModifierBase.h"

namespace qmcplusplus
{
DriftModifierBase* createDriftModifier(xmlNodePtr cur);
} // namespace qmcplusplus

#endif
