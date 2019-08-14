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

#include "QMCDrivers/VMC/VMCDriverInput.h"

namespace qmcplusplus
{
VMCDriverInput::VMCDriverInput(int walkers_per_rank, const std::string& use_drift)
    : requested_walkers_per_rank_(walkers_per_rank), use_drift_(use_drift)
{}

void VMCDriverInput::readXML(xmlNodePtr& node)
{
  ParameterSet parameter_set_;
  parameter_set_.add(requested_walkers_per_rank_, "walkers", "int");
  parameter_set_.add(use_drift_, "usedrift", "string");
  parameter_set_.add(use_drift_, "use_drift", "string");

  parameter_set_.put(node);
}

} // namespace qmcplusplus
