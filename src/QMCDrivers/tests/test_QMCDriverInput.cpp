//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
#include <catch2/catch_test_macros.hpp>
#include "Utilities/for_testing/Catch2Approx.h"

#include "QMCDrivers/QMCDriverInput.h"
#include "EstimatorInputDelegates.h"
#include "QMCDrivers/tests/ValidQMCInputSections.h"
#include "OhmmsData/Libxml2Doc.h"

namespace qmcplusplus
{
TEST_CASE("QMCDriverInput Instantiation", "[drivers]") { QMCDriverInput driver_input; }

TEST_CASE("QMCDriverInput readXML", "[drivers]")
{
  auto xml_test = [](const char* driver_xml) {
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(driver_xml));
    xmlNodePtr node = doc.getRoot();
    QMCDriverInput qmcdriver_input;
    qmcdriver_input.readXML(node);
    REQUIRE(qmcdriver_input.get_qmc_method().size() > 0);
  };

  std::for_each(testing::valid_vmc_input_sections.begin() + testing::valid_vmc_input_vmc_batch_index,
                testing::valid_vmc_input_sections.end(), xml_test);

  std::for_each(testing::valid_dmc_input_sections.begin() + testing::valid_dmc_input_dmc_batch_index,
                testing::valid_dmc_input_sections.end(), xml_test);
}

TEST_CASE("QMCDriverInput retired check-properties parameters", "[drivers]")
{
  const int default_blocks_between_recompute =
      std::is_same<QMCTraits::RealType, QMCTraits::FullPrecRealType>::value ? 10 : 1;

  SECTION("omitted")
  {
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(R"(<qmc method="vmc" move="pbyp"/>)"));
    QMCDriverInput qmcdriver_input;
    qmcdriver_input.readXML(doc.getRoot());
    CHECK(qmcdriver_input.get_blocks_between_recompute() == default_blocks_between_recompute);
    CHECK(qmcdriver_input.get_debug_checks() == DriverDebugChecks::ALL_OFF);
  }

  for (const std::string& alias : {"checkProperties", "checkproperties", "check_properties"})
  {
    CAPTURE(alias);
    const std::string input =
        "<qmc method=\"vmc\" move=\"pbyp\">"
        "<parameter name=\"" +
        alias +
        "\">not-an-integer</parameter>"
        "<parameter name=\"blocks_between_recompute\">3</parameter>"
        "<parameter name=\"debug_checks\">checkGL_after_moves</parameter>"
        "</qmc>";
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(input));
    QMCDriverInput qmcdriver_input;
    qmcdriver_input.readXML(doc.getRoot());
    CHECK(qmcdriver_input.get_blocks_between_recompute() == 3);
    CHECK(qmcdriver_input.get_debug_checks() == DriverDebugChecks::CHECKGL_AFTER_MOVES);
  }
}

TEST_CASE("QMCDriverInput retired modern state", "[drivers]")
{
  Libxml2Document doc;
  REQUIRE(doc.parseFromString(R"(<qmc method="vmc" move="pbyp">
      <parameter name="current">not-an-integer</parameter>
      <parameter name="maxDisplSq">not-a-real</parameter>
      <parameter name="blocks_between_recompute">3</parameter>
      <random seed="retired-qmc-child"/>
    </qmc>)"));

  QMCDriverInput qmcdriver_input;
  qmcdriver_input.readXML(doc.getRoot());
  CHECK(qmcdriver_input.get_blocks_between_recompute() == 3);
  CHECK(qmcdriver_input.get_qmc_method() == "vmc");
}

TEST_CASE("QMCDriverInput retired configuration dump controls", "[drivers]")
{
  Libxml2Document doc;
  REQUIRE(doc.parseFromString(R"(<qmc method="vmc" move="pbyp" checkpoint="3">
      <parameter name="recordconfigs">not-an-integer</parameter>
      <parameter name="record_configs">also-not-an-integer</parameter>
      <record stride="not-an-integer" period="also-not-an-integer"/>
      <checkpoint stride="not-an-integer" period="4"/>
      <dumpconfig stride="not-an-integer" period="also-not-an-integer"/>
    </qmc>)"));

  QMCDriverInput qmcdriver_input;
  qmcdriver_input.readXML(doc.getRoot());
  CHECK(qmcdriver_input.get_check_point_period() == 4);
  CHECK(qmcdriver_input.get_dump_config());
}
} // namespace qmcplusplus
