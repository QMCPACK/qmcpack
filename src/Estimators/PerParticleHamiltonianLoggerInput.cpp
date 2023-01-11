//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File crreated by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "PerParticleHamiltonianLoggerInput.h"
#include "EstimatorInput.h"

namespace qmcplusplus {
  PerParticleHamiltonianLoggerInput::PerParticleHamiltonianLoggerInput(xmlNodePtr cur)
  {
    input_section_.readXML(cur);
    auto setIfInInput = LAMBDA_setIfInInput;
    setIfInInput(to_stdout_, "to_stdout");
    setIfInInput(name_, "name");
  }
}
