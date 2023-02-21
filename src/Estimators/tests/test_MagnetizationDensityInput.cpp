
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
  using PosType = QMCTraits::PosType;
  using namespace testing::magdensity;

  Lattice lattice = testing::makeTestLattice();

  for (auto input_xml : valid_mag_density_input_sections)
  {
    Libxml2Document doc;
    bool okay = doc.parseFromString(input_xml);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    MagnetizationDensityInput magdens(node);
    MagnetizationDensityInput::DerivedParameters dev_par = magdens.calculateDerivedParameters(lattice);
  }

  for (auto input_xml : testing::invalid_mag_density_input_sections)
  {
    Libxml2Document doc;
    bool okay = doc.parseFromString(input_xml);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();

    CHECK_THROWS(MagnetizationDensityInput(node));
  }

  //Lastly, going to check the state variables are consistent with the input.
  //We'll take the magdensity::Inputs::valid_magdensity_input test case for
  //thorough checking.
  {
    Libxml2Document doc;
    auto input_xml = valid_mag_density_input_sections[Inputs::valid_magdensity_input];
    bool okay      = doc.parseFromString(input_xml);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    MagnetizationDensityInput magdens(node);
    int nsamples   = magdens.get_nsamples();
    int integrator = int(magdens.get_integrator());
    PosType grid   = magdens.get_grid();
    PosType dr     = magdens.get_dr();
    PosType corner = magdens.get_corner();


    app_log() << "NSAMPLES = " << nsamples << std::endl;
    app_log() << "INTEGRATOR = " << integrator << std::endl;
    app_log() << "CORNER = " << corner << std::endl;
    app_log() << "GRID = " << grid << std::endl;
    app_log() << "DR   = " << dr << std::endl;

    CHECK(nsamples == 64);
    CHECK(integrator == 0);
    CHECK(grid[0] == Approx(4));
    CHECK(grid[1] == Approx(3));
    CHECK(grid[2] == Approx(2));
    CHECK(corner[0] == Approx(0.0));
    CHECK(corner[1] == Approx(0.0));
    CHECK(corner[2] == Approx(0.0));
    CHECK(dr[0] == Approx(0.1));
    CHECK(dr[1] == Approx(0.1));
    CHECK(dr[2] == Approx(0.1));
  }
}

} // namespace qmcplusplus
