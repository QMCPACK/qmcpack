//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "DMCDriverInput.h"

namespace qmcplusplus
{
DMCDriverInput::DMCDriverInput(int walkers_per_rank) {}

void DMCDriverInput::readXML(xmlNodePtr node)
{
  ParameterSet parameter_set_;
  std::string reconfig_str;
  parameter_set_.add(reconfig_str, "reconfiguration");
  if (!reconfig_str.empty() && reconfig_str != "no" && reconfig_str != "runwhileincorrect")
    throw std::runtime_error(
        "Reconfiguration is currently broken and gives incorrect results. Set reconfiguration=\"no\" or remove the "
        "reconfiguration option from the DMC input section. To run performance tests, please set reconfiguration to "
        "\"runwhileincorrect\" instead of \"yes\" to restore consistent behaviour.");
  reconfiguration_ = (reconfig_str == "runwhileincorrect");
  parameter_set_.add(NonLocalMove, "nonlocalmove");
  parameter_set_.add(NonLocalMove, "nonlocalmoves");
  parameter_set_.add(max_age_, "MaxAge");

  // from DMC.cpp put(xmlNodePtr)
  parameter_set_.add(branch_interval_, "branchInterval");
  parameter_set_.add(branch_interval_, "branchinterval");
  parameter_set_.add(branch_interval_, "substeps");
  parameter_set_.add(branch_interval_, "subStep");
  parameter_set_.add(branch_interval_, "sub_stepd");

  //from NonLocalTOperator.cpp
  parameter_set_.add(alpha_, "alpha");
  parameter_set_.add(gamma_, "gamma");

  parameter_set_.add(reserve_, "reserve");

  parameter_set_.put(node);

  // TODO: similar check for alpha and gamma
  if (max_age_ < 0)
    throw std::runtime_error("Illegal input for MaxAge in DMC input section");
  if (branch_interval_ < 0)
    throw std::runtime_error("Illegal input for branchInterval or substeps in DMC input section");

  if (reserve_ < 1.0)
    throw std::runtime_error("You can only reserve walkers above the target walker count");
}

std::ostream& operator<<(std::ostream& o_stream, const DMCDriverInput& dmci) { return o_stream; }

} // namespace qmcplusplus
