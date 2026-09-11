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
#include <catch2/catch_test_macros.hpp>
#include "Utilities/for_testing/Catch2Approx.h"

#include "OneBodyDensityMatricesInput.h"
#include "ValidOneBodyDensityMatricesInput.h"
#include "InvalidOneBodyDensityMatricesInput.h"
#include "EstimatorTesting.h"
#include "ParticleSet.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Message/UniformCommunicateError.h"

#include <iostream>
#include <string_view>

namespace qmcplusplus
{
TEST_CASE("OneBodyDensityMatricesInput::from_xml", "[estimators]")
{
  using Input = testing::ValidOneBodyDensityMatricesInput;
  Input valid_input;
  for (auto input_xml : valid_input)
  {
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(input_xml));
    xmlNodePtr node = doc.getRoot();
    OneBodyDensityMatricesInput obdmi(node);
  }

  using invalid_input = testing::InvalidOneBodyDensityMatricesInput;
  for (auto input_xml : invalid_input::xml)
  {
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(input_xml));
    xmlNodePtr node = doc.getRoot();

    CHECK_THROWS_AS(OneBodyDensityMatricesInput(node), UniformCommunicateError);
  }
}

TEST_CASE("OneBodyDensityMatricesInput::copy_construction", "[estimators]")
{
  using Input = testing::ValidOneBodyDensityMatricesInput;
  Libxml2Document doc;
  REQUIRE(doc.parseFromString(Input::getXml(Input::valid::SCALE)));
  xmlNodePtr node = doc.getRoot();
  OneBodyDensityMatricesInput obdmi(node);
  static_assert(std::is_copy_constructible_v<OneBodyDensityMatricesInput>);
}

TEST_CASE("OneBodyDensityMatricesInput::volume_normed", "[estimators]")
{
  constexpr std::string_view input_xml = R"XML(
<estimator type="OneBodyDensityMatrices" name="OneBodyDensityMatrices">
  <parameter name="basis">spo_ud</parameter>
  <parameter name="integrator">uniform</parameter>
  <parameter name="volume_normed">no</parameter>
</estimator>
)XML";

  Libxml2Document doc;
  REQUIRE(doc.parseFromString(input_xml));
  OneBodyDensityMatricesInput obdmi(doc.getRoot());
  CHECK_FALSE(obdmi.get_volume_normalized());
}

} // namespace qmcplusplus
