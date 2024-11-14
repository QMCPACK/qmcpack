//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "StructureFactorInput.h"
#include "ValidStructureFactorInput.h"
#include "OhmmsData/Libxml2Doc.h"
#include <iostream>

namespace qmcplusplus
{

TEST_CASE("StructureFactorInput::parseXML::valid", "[estimators]")
{
  using input = qmcplusplus::testing::ValidStructureFactorInput;
  int test_num = 0;
  for (auto input_xml : input::xml)
  {
    std::cout << "input number: " << test_num++ << '\n'; 
    Libxml2Document doc;
    bool okay       = doc.parseFromString(input_xml);
    xmlNodePtr node = doc.getRoot();
    StructureFactorInput sfi(node);
  }
}

TEST_CASE("StructureFactorInput::parseXML::invalid", "[estimators]")
{
  using input = qmcplusplus::testing::InvalidStructureFactorInput;
  for (auto input_xml : input::xml)
  {
    Libxml2Document doc;
    bool okay                           = doc.parseFromString(input_xml);
    xmlNodePtr node                     = doc.getRoot();
    auto constructBadStructureFactorInput = [](xmlNodePtr cur) { StructureFactorInput sfi(cur); };
    CHECK_THROWS_AS(constructBadStructureFactorInput(node), UniformCommunicateError);
  }
}

}
