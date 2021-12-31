//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: EstimatorManagerNew.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "EstimatorManagerInput.h"
#include "MomentumDistributionInput.h"
#include "OneBodyDensityMatricesInput.h"
#include "SpinDensityInput.h"

namespace qmcplusplus
{

EstimatorManagerInput::EstimatorManagerInput(xmlNodePtr cur)
{
  readXML(cur);
}

void EstimatorManagerInput::readXML(xmlNodePtr cur)
{
  xmlNodePtr child = cur->xmlChildrenNode;
  while(child != NULL)
  {
  }
}


}


