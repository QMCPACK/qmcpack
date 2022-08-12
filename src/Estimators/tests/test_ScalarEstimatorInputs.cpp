//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "ScalarEstimatorInputs.h"
#include "ValidScalarEstimatorInput.h"
#include "ModernStringUtils.hpp"

namespace qmcplusplus
{

TEST_CASE("Scalar Estimator Input", "[estimators]")
{
  using namespace testing;

  for (auto input : valid_scalar_estimator_input_sections)
  {
    Libxml2Document doc;
    bool okay = doc.parseFromString(input);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    std::string atype(lowerCase(getXMLAttributeValue(node, "type")));
    std::string aname(lowerCase(getXMLAttributeValue(node, "name")));
    // Since legacy inconsistently used name instead of type attribute to specify scalar estimator type
    if (atype.empty() && !aname.empty())
      atype = aname;
    if (aname.empty() && !atype.empty())
      aname = atype;
    if (atype == "localenergy" || atype == "elocal")
    {
      LocalEnergyInput lei(node);
      CHECK(lowerCase(lei.get_type()) == "localenergy");
    }
    else if (atype == "cslocalenergy")
    {
      CSLocalEnergyInput cslei(node);
      CHECK(lowerCase(cslei.get_type()) == "cslocalenergy");
    }
    else if (atype == "rmc")
    {
      RMCLocalEnergyInput rmclei(node);
      CHECK(lowerCase(rmclei.get_type()) == "rmclocalenergy");
    }
    else
    {
      FAIL("Unhandled scalar estimator  " << atype);
    }
  }
}


} // namespace qmcplusplus
