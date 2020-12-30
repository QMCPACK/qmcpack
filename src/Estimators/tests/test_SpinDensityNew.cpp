//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "SpinDensityInput.h"
#include "ValidSpinDensityInput.h"
#include "SpinDensityNew.h"
#include "Utilities/SpeciesSet.h"

#include "OhmmsData/Libxml2Doc.h"

#include <stdio.h>
#include <sstream>

namespace qmcplusplus
{
TEST_CASE("SpinDensityNew", "[estimators]")
{
  Libxml2Document doc;
  bool okay = doc.parseFromString(testing::valid_spin_density_input_sections[0]);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();
  SpinDensityInput sdi;
  sdi.readXML(node);
  SpeciesSet species_set;
  int ispecies = species_set.addSpecies("C");
  int iattribute = species_set.addAttribute("membersize");
  species_set(iattribute, ispecies) = 2;
  SpinDensityNew(sdi, species_set);
}

TEST_CASE("SpinDensityNew::accumulate", "[estimators]")
{
  Libxml2Document doc;
  bool okay = doc.parseFromString(testing::valid_spin_density_input_sections[0]);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();
  SpinDensityInput sdi;
  sdi.readXML(node);
  SpeciesSet species_set;
  int ispecies = species_set.addSpecies("C");
  int iattribute = species_set.addAttribute("membersize");
  species_set(iattribute, ispecies) = 2;
  SpinDensityNew(sdi, species_set);
}

TEST_CASE("SpinDensityNew::collect", "[estimators]") {}
}
