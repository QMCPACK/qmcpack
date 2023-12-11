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

#include "QMCDrivers/GreenFunctionModifiers/DriftModifierBase.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{
/// create DriftModifier
DriftModifierBase* createDriftModifier(xmlNodePtr cur, const Communicate* myComm);
DriftModifierBase* createDriftModifier(const std::string& drift_modifier_str, QMCTraits::RealType unr_a);
} // namespace qmcplusplus

#endif
