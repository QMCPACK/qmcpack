//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, mdewing@anl.gov Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"


//#include "Utilities/RandomGenerator.h"
#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsApp/RandomNumberControl.h"


#include <stdio.h>
#include <string>



namespace qmcplusplus
{

TEST_CASE("RandomNumberControl make_seeds", "[ohmmsapp]")
{

  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  RandomNumberControl::make_seeds();

  REQUIRE(RandomNumberControl::Children.size() > 0);
}

TEST_CASE("RandomNumberControl no random in xml", "[ohmmsapp]")
{
  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  const char *xml_input="<tmp></tmp>";

  Libxml2Document doc;
  bool okay = doc.parseFromString(xml_input);
  REQUIRE(okay);

  RandomNumberControl rnc;

  xmlXPathContextPtr context = doc.getXPathContext();
  rnc.initialize(context);
}

TEST_CASE("RandomNumberControl random in xml", "[ohmmsapp]")
{
  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  const char *xml_input="<tmp><random seed='0'></random></tmp>";

  Libxml2Document doc;
  bool okay = doc.parseFromString(xml_input);
  REQUIRE(okay);

  RandomNumberControl rnc;

  xmlXPathContextPtr context = doc.getXPathContext();
  rnc.initialize(context);


  rnc.write("rng_out",c);

  RandomNumberControl rnc2;
  rnc2.read("rng_out",c);
  // not sure what to test here - for now make sure it doesn't crash.
}
}

