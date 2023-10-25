//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include <stdio.h>
#include <sstream>

#include "ReferencePointsInput.h"
#include "ValidReferencePointsInput.h"
#include "OhmmsData/Libxml2Doc.h"

namespace qmcplusplus
{

TEST_CASE("ReferencePointsInput::parseXML::valid", "[estimators]")
{
  using Input = testing::ValidReferencePointsInputs;
  for (auto input_xml : Input::xml)
  {
    Libxml2Document doc;
    bool okay       = doc.parseFromString(input_xml);
    xmlNodePtr node = doc.getRoot();

    // Will throw if input is invalid.
    ReferencePointsInput rpi(node);
  }
}

TEST_CASE("ReferencePointsInput::parseXML::invalid", "[estimators]")
{
  using Input = testing::InvalidReferencePointsInputs;
  for (auto input_xml : Input::xml)
  {
    Libxml2Document doc;
    bool okay = doc.parseFromString(input_xml);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();

    auto constructBadRefPoints = [](xmlNodePtr cur) { ReferencePointsInput rpi(cur); };
    CHECK_THROWS_AS(constructBadRefPoints(node), UniformCommunicateError);
  }
}

TEST_CASE("ReferencePointsInput::makeReferencePointsInput", "[estimators]")
{
  using Input = testing::ValidReferencePointsInputs;
  std::string value_label;
  Libxml2Document doc;
  bool okay       = doc.parseFromString(Input::xml[0]);
  xmlNodePtr node = doc.getRoot();
  makeReferencePointsInput(node, value_label);
  CHECK(value_label == "referencepoints");
}
} // namespace qmcplusplus
