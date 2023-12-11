//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: MomentumDistribution.cpp
//////////////////////////////////////////////////////////////////////////////////////
#include "MomentumDistributionInput.h"

namespace qmcplusplus
{

MomentumDistributionInput::MomentumDistributionInput(xmlNodePtr cur)
{
  // This results in checkParticularValidity being called on MomentumDistributionInput
  input_section_.readXML(cur);

  auto setIfInInput = [&](auto& var, const std::string& tag) -> bool { return input_section_.setIfInInput(var, tag); };
  setIfInInput(name_, "name");
  setIfInInput(type_, "type");
  setIfInInput(samples_, "samples");
  setIfInInput(kmax_, "kmax");
  setIfInInput(kmax0_, "kmax0");
  setIfInInput(kmax1_, "kmax1");
  setIfInInput(kmax2_, "kmax2");
}


} // namespace qmcplusplus

