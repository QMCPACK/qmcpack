
//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories 
//
// File created by:  Raymond Clay, rclay@sandia.gov, Sandia National Laboratories 
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "MagnetizationDensityInput.h"
#include "ValidMagnetizationDensityInput.h"
#include "InvalidMagnetizationDensityInput.h"
#include "EstimatorTesting.h"
#include "ParticleSet.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Message/UniformCommunicateError.h"

#include <iostream>

namespace qmcplusplus
{
TEST_CASE("MagnetizationDensityInput::from_xml", "[estimators]")
{
  using POLT    = PtclOnLatticeTraits;
  using Lattice = POLT::ParticleLayout;
  using namespace testing::magdensity;
  MagnetizationDensityInput magdens_in;

  for (auto input_xml : valid_mag_density_input_sections)
  {
    Libxml2Document doc;
    bool okay = doc.parseFromString(input_xml);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    MagnetizationDensityInput obdmi(node);
  }

  for (auto input_xml : testing::invalid_mag_density_input_sections)
  {
    Libxml2Document doc;
    bool okay = doc.parseFromString(input_xml);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();

  //  CHECK_THROWS_AS(MagnetizationDensityInput(node), UniformCommunicateError);
  }
  
}

} // namespace qmcplusplus
