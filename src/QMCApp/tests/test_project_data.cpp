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


#include "Message/catch_mpi_main.hpp"


#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsApp/ProjectData.h"


#include <stdio.h>
#include <string>
#include <sstream>


using std::string;

namespace qmcplusplus
{

TEST_CASE("ProjectData", "[ohmmsapp]")
{

  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;


  ProjectData proj1;
  // If no name given, it gets set to time and date
  //   and the title is set equal to the name
  REQUIRE(std::string(proj1.m_title).size() > 0);
  REQUIRE(std::string(proj1.m_title) == proj1.getName());


  ProjectData proj2("test");
  REQUIRE(proj2.m_series == 0);
  proj2.advance();
  REQUIRE(proj2.m_series == 1);

  REQUIRE(proj2.m_title == std::string("asample"));
  REQUIRE(proj2.getName() == std::string("test"));

  proj2.setCommunicator(c);
  std::stringstream o2;
  proj2.get(o2);
}

TEST_CASE("ProjectData::put no series", "[ohmmsapp]")
{

  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  ProjectData proj("test");

  const char *xml_input = "<project id='test1'></project>";
  Libxml2Document doc;
  bool okay = doc.parseFromString(xml_input);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  proj.put(root);
  REQUIRE(proj.m_series == 0);
}

TEST_CASE("ProjectData::put with series", "[ohmmsapp]")
{

  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  ProjectData proj("test");

  const char *xml_input = "<project id='test1' series='1'></project>";
  Libxml2Document doc;
  bool okay = doc.parseFromString(xml_input);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  proj.put(root);
  REQUIRE(proj.m_series == 1);

  // host and date nodes get added for output to the .cont.xml file

}

}

