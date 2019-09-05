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
VMCDriverInput::VMCDriverInput(int walkers_per_rank, bool use_drift)
    : walkers_per_rank_(walkers_per_rank), use_drift_(use_drift)
{}

void VMCDriverInput::readXML(xmlNodePtr node)
{
  ParameterSet parameter_set_;
  parameter_set_.add(walkers_per_rank_, "walkers", "int");
  std::string use_drift;
  parameter_set_.add(use_drift, "usedrift", "string");
  parameter_set_.add(use_drift, "use_drift", "string");
  parameter_set_.add(samples_, "samples", "int");
  parameter_set_.add(samples_per_thread_, "samplesperthread", "int");
  parameter_set_.add(steps_between_samples_, "stepsbetweensamples", "int");
  parameter_set_.put(node);
  if(use_drift == "yes")
    use_drift_ = true;
}

std::ostream& operator<<(std::ostream& o_stream, const VMCDriverInput& vmci) { return o_stream; }

} // namespace qmcplusplus
