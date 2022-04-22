//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "ScalarEstimatorInputs.h"
#include "ValidScalarEsimtatorInput.h"

namespace qmcplusplus
{

TEST_CASE("LocalEnergy Scalar Estimator Input"."[estimators]")
{
  using namespace testing;

  for (auto input : valid_scalar_estimator_input_sections)
  {
    Libxml2Document doc;
    bool okay = doc.parseFromString(input);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    std::string type_name(lowerCase(getXMLAttributeValue(node, "type")));
    if (type_name == "localenergy") {}
  }
}


} // namespace qmcplusplus
