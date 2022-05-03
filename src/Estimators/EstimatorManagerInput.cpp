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
  std::string cur_name{lowerCase(castXMLCharToChar(cur->name))};
  xmlNodePtr child;
  if (cur_name == "estimators")
    child = cur->xmlChildrenNode;
  else
    child = cur;
  while (child != NULL)
  {
    std::string cname{lowerCase(castXMLCharToChar(child->name))};
    if (cname == "estimator")
    {
      std::string atype(lowerCase(getXMLAttributeValue(child, "type")));
      std::string aname(lowerCase(getXMLAttributeValue(child, "name")));
      if (atype.empty() && !aname.empty())
        atype = aname;
      if (aname.empty() && !atype.empty())
        aname = atype;
      if (atype == "localenergy" || atype == "elocal")
        appendScalarEstimatorInput<LocalEnergyInput>(child);
      else if (atype == "cslocalenergy")
      {
        appendScalarEstimatorInput<CSLocalEnergyInput>(child);
        app_warning() << "CSLocalEnergyEstimator support is at best experimental with batch drivers";
      }
      else if (atype == "rmc")
      {
        appendScalarEstimatorInput<RMCLocalEnergyInput>(child);
        app_warning() << "RMCLocalEnergyEstimator support is at best experimental with batch drivers";
      }
      else if (atype == "onebodydensitymatrices")
        appendEstimatorInput<OneBodyDensityMatricesInput>(child);
      else if (atype == "spindensity")
        appendEstimatorInput<SpinDensityInput>(child);
      else if (atype == "momentumdistribution")
        appendEstimatorInput<MomentumDistributionInput>(child);
      else
        throw UniformCommunicateError(error_tag + "unparsable <estimator> node, name: " + aname + " type: " + atype +
                                      " in Estimators input.");
    }
    else if (cname != "text")
    {
      std::string atype(lowerCase(getXMLAttributeValue(child, "type")));
      std::string aname(lowerCase(getXMLAttributeValue(child, "name")));
      throw UniformCommunicateError(error_tag + "<Estimators> can only contain <Estimator> nodes");
    }

    if (cur_name == "estimators")
      child = child->next;
    else
    {
      app_summary() << "<estimator> nodes not contained in <estimators></estimators> is a deprecated input xml idiom";
      break;
    }
  }
}


} // namespace qmcplusplus
