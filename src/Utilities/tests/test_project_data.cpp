//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
#include <catch2/catch_test_macros.hpp>


#include "OhmmsData/Libxml2Doc.h"
#include "ProjectData.h"
#include "Utilities/qmc_common.h"


#include <stdio.h>
#include <string>
#include <sstream>


using std::string;

namespace qmcplusplus
{
TEST_CASE("ProjectData", "[ohmmsapp]")
{
  ProjectData proj1;
  // If no name given, it gets set to time and date
  //   and the title is set equal to the name
  REQUIRE(proj1.getTitle().size() > 0);

  ProjectData proj2("asample");
  REQUIRE(proj2.getSeriesIndex() == 0);
  proj2.advance();
  REQUIRE(proj2.getSeriesIndex() == 1);
  REQUIRE(proj2.getTitle() == std::string("asample"));

  std::stringstream o2;
  proj2.get(o2);

  qmc_common.setupFilenameInfixForEnsembleJobs(1, 4);
  ProjectData proj_multi("my_group");
  REQUIRE(proj_multi.currentMainRoot() == "my_group.g001.s000");

  qmc_common.setupFilenameInfixForEnsembleJobs(0, 1);
  ProjectData proj_single("single_group");
  REQUIRE(proj_single.currentMainRoot() == "single_group.s000");
}

TEST_CASE("ProjectData::put no series", "[ohmmsapp]")
{
  ProjectData proj;

  const char* xml_input = R"(<project id="test1"></project>)";
  Libxml2Document doc;
  REQUIRE(doc.parseFromString(xml_input));

  xmlNodePtr root = doc.getRoot();

  proj.put(root);
  REQUIRE(proj.getSeriesIndex() == 0);
  REQUIRE(proj.getTitle() == std::string("test1"));
}

TEST_CASE("ProjectData::put with series", "[ohmmsapp]")
{
  ProjectData proj;

  const char* xml_input = R"(<project id="test1" series="1"></project>)";
  Libxml2Document doc;
  REQUIRE(doc.parseFromString(xml_input));

  xmlNodePtr root = doc.getRoot();

  proj.put(root);
  REQUIRE(proj.getSeriesIndex() == 1);
  REQUIRE(proj.getTitle() == std::string("test1"));

  // host and date nodes get added for output to the .cont.xml file
}

TEST_CASE("ProjectData::TestDriverVersion", "[ohmmsapp]")
{
  using DV = ProjectData::DriverVersion;
  SECTION("driver version batch")
  {
    ProjectData proj_batch(DV::BATCH);
    REQUIRE(proj_batch.getDriverVersion() == DV::BATCH);

    ProjectData proj;
    REQUIRE(proj.getDriverVersion() == DV::BATCH);

    const char* xml_input = R"(
      <project id="test1" series="1">
        <parameter name='driver_version'>
          batch
        </parameter>
      </project>
      )";
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(xml_input));

    xmlNodePtr root = doc.getRoot();

    proj.put(root);
    REQUIRE(proj.getSeriesIndex() == 1);
    REQUIRE(proj.getDriverVersion() == DV::BATCH);
  }
  SECTION("driver version legacy")
  {
    ProjectData proj_legacy(DV::LEGACY);
    REQUIRE(proj_legacy.getDriverVersion() == DV::LEGACY);

    ProjectData proj;
    const char* xml_input = R"(
      <project id="test1" series="1">
        <parameter name='driver_version'>
          legacy
        </parameter>
      </project>
      )";
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(xml_input));

    xmlNodePtr root = doc.getRoot();

    proj.put(root);
    REQUIRE(proj.getSeriesIndex() == 1);
    REQUIRE(proj.getDriverVersion() == DV::LEGACY);
  }
  SECTION("driver version bad value")
  {
    ProjectData proj;

    const char* xml_input = R"(
      <project id=" test1 " series=" 1 ">
        <parameter name = 'driver_version' >
          linear
        </parameter>
      </project>)";
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(xml_input));

    xmlNodePtr root = doc.getRoot();

    REQUIRE_THROWS(proj.put(root));
  }

  // host and date nodes get added for output to the .cont.xml file
}

TEST_CASE("ProjectData::TestIsComplex", "[ohmmsapp]")
{
  ProjectData proj;
#ifdef QMC_COMPLEX
  REQUIRE(proj.isComplex());
#else
  REQUIRE(!proj.isComplex());
#endif
}

} // namespace qmcplusplus
