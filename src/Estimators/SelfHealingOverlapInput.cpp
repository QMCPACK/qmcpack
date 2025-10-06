//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "SelfHealingOverlapInput.h"
#include "EstimatorInput.h"

namespace qmcplusplus
{
SelfHealingOverlapInput::SelfHealingOverlapInput(xmlNodePtr cur)
{
  input_section_.readXML(cur);
  auto setIfInInput = LAMBDA_setIfInInput;
  setIfInInput(name_, "name");
  setIfInInput(type_, "type");
}

} // namespace qmcplusplus
