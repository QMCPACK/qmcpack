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

#include "VMCDriverInput.h"

namespace qmcplusplus
{
VMCDriverInput::VMCDriverInput(bool use_drift) : use_drift_(use_drift) {}

void VMCDriverInput::readXML(xmlNodePtr node)
{
  ParameterSet parameter_set_;
  std::string use_drift;
  parameter_set_.add(use_drift, "usedrift", {"yes", "no"});
  parameter_set_.add(use_drift, "use_drift", {"yes", "no"});
  parameter_set_.add(samples_, "samples");
  parameter_set_.put(node);

  use_drift_ = use_drift == "yes";
  if (use_drift_)
    app_log() << "  Random walking with drift" << std::endl;
  else
    app_log() << "  Random walking without drift" << std::endl;
}

std::ostream& operator<<(std::ostream& o_stream, const VMCDriverInput& vmci) { return o_stream; }

} // namespace qmcplusplus
