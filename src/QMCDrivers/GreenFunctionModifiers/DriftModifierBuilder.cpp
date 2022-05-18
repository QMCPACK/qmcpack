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


#include "DriftModifierBuilder.h"
#include "QMCDrivers/GreenFunctionModifiers/DriftModifierUNR.h"
#include "OhmmsData/ParameterSet.h"

namespace qmcplusplus
{
DriftModifierBase* createDriftModifier(xmlNodePtr cur, const Communicate* myComm)
{
  std::string ModifierName("UNR");
  ParameterSet m_param;
  m_param.add(ModifierName, "drift_modifier");
  m_param.put(cur);
  if (ModifierName != "UNR")
    myComm->barrier_and_abort("createDriftModifier unknown drift_modifier " + ModifierName);
  DriftModifierBase* DriftModifier = new DriftModifierUNR;
  return DriftModifier;
}

DriftModifierBase* createDriftModifier(const std::string& drift_modifier_str, QMCTraits::RealType unr_a)
{
  std::string dm_str(lowerCase(drift_modifier_str));
  if (dm_str != "unr")
    throw std::runtime_error("createDriftModifier unknown drift_modifier ");
  DriftModifierBase* DriftModifier = new DriftModifierUNR(unr_a);
  return DriftModifier;
}

} // namespace qmcplusplus
