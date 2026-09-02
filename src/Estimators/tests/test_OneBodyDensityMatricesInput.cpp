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
#include <catch2/matchers/catch_matchers_string.hpp>
#include "Utilities/for_testing/Catch2Approx.h"

#include "OneBodyDensityMatricesInput.h"
#include "ValidOneBodyDensityMatricesInput.h"
#include "InvalidOneBodyDensityMatricesInput.h"
#include "EstimatorTesting.h"
#include "ParticleSet.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Message/UniformCommunicateError.h"

#include <array>
#include <iostream>
#include <string>
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

TEST_CASE("OneBodyDensityMatricesInput rejects legacy-only diagnostics", "[estimators]")
{
  constexpr std::array<std::string_view, 3> legacy_only_diagnostics{
      "energy_matrix", "check_overlap", "check_derivatives"};

  for (const std::string_view diagnostic : legacy_only_diagnostics)
  {
    DYNAMIC_SECTION("reject " << diagnostic)
    {
      const std::string input_xml =
          "<estimator type=\"OneBodyDensityMatrices\" name=\"OneBodyDensityMatrices\">"
          "<parameter name=\"basis\">spo_ud</parameter>"
          "<parameter name=\"integrator\">uniform</parameter>"
          "<parameter name=\"" +
          std::string(diagnostic) +
          "\">yes</parameter>"
          "</estimator>";

      Libxml2Document doc;
      REQUIRE(doc.parseFromString(input_xml));
      CHECK_THROWS_WITH(OneBodyDensityMatricesInput(doc.getRoot()),
                        Catch::Matchers::ContainsSubstring("name " + std::string(diagnostic) + " is not a parameter"));
    }
  }
}

} // namespace qmcplusplus
