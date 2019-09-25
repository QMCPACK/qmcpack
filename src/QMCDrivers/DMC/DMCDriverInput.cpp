//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCDrivers/DMC/DMCDriverInput.h"

namespace qmcplusplus
{
DMCDriverInput::DMCDriverInput(int walkers_per_rank) {}

void DMCDriverInput::readXML(xmlNodePtr node)
{
  ParameterSet parameter_set_;
  std::string reconfig_str;
  parameter_set_.add(reconfig_str, "reconfiguration", "string");
  reconfiguration = (reconfig_str == "yes");
  parameter_set_.add(NonLocalMove, "nonlocalmove", "string");
  parameter_set_.add(NonLocalMove, "nonlocalmoves", "string");
  parameter_set_.add(mover_MaxAge, "MaxAge", "double");

  // from DMC.cpp put(xmlNodePtr)
  parameter_set_.add(BranchInterval, "branchInterval", "string");
  parameter_set_.add(BranchInterval, "branchinterval", "string");
  parameter_set_.add(BranchInterval, "substeps", "int");
  parameter_set_.add(BranchInterval, "subStep", "int");
  parameter_set_.add(BranchInterval, "sub_stepd", "int");

  parameter_set_.put(node);
}

std::ostream& operator<<(std::ostream& o_stream, const DMCDriverInput& dmci) { return o_stream; }

} // namespace qmcplusplus
