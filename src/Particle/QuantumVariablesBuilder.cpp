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


#include "QuantumVariablesBuilder.h"
#include "Particle/RealSpacePositions.h"

namespace qmcplusplus
{
/** create QuantumVariables based on kind
 */
std::unique_ptr<QuantumVariables> createQuantumVariables(const QuantumVariableKind kind)
{
  if (kind == QuantumVariableKind::QV_POS)
    return std::make_unique<RealSpacePositions>();
  else
    APP_ABORT("QuantumVariablesBuilder::createQuantumVariables unknown QuantumVariableKind");
  // dummy return
  return std::unique_ptr<RealSpacePositions>();
}
} // namespace qmcplusplus
