//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: EstimatorManagerNew.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "EstimatorManagerInput.h"
#include "ScalarEstimatorInputs.h"
#include "MomentumDistributionInput.h"
#include "OneBodyDensityMatricesInput.h"
#include "SpinDensityInput.h"

namespace qmcplusplus
{
EstimatorManagerInput::EstimatorManagerInput(xmlNodePtr cur) { readXML(cur); }

void EstimatorManagerInput::readXML(xmlNodePtr cur)
{
  const std::string error_tag{"EstimatorManager input:"};
  xmlNodePtr child = cur->xmlChildrenNode;
  while (child != NULL)
  {
    std::string cname{lowerCase(castXMLCharToChar(child->name))};
    if (cname == "estimator")
    {
      std::string atype(lowerCase(getXMLAttributeValue(child, "type")));
      if (atype == "localenergy")
        appendScalarEstimatorInput<LocalEnergyInput>(child);
      else if (atype == "cslocalenergy")
        appendScalarEstimatorInput<CSLocalEnergyInput>(child);
      else if (atype == "onebodydensitymatrices")
        appendEstimatorInput<OneBodyDensityMatricesInput>(child);
      else if (atype == "spindensity")
        appendEstimatorInput<SpinDensityInput>(child);
      else if (atype == "momentumdistribution")
        appendEstimatorInput<MomentumDistributionInput>(child);
      else
        throw UniformCommunicateError(error_tag + "unparsable child node, name: " + cname + " type: " + atype +
                                      " in Estimators input.");
    }
    else
      throw UniformCommunicateError(error_tag + "<Estimators> can only contain <Estimator> nodes");
    child = child->next;
  }
}


} // namespace qmcplusplus
