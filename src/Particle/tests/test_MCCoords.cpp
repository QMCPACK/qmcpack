//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "MCCoords.hpp"
#include "Utilities/StlPrettyPrint.hpp"
#include "Utilities/StdRandom.h"
#include "ParticleBase/RandomSeqGenerator.h"

namespace qmcplusplus
{

TEST_CASE("MCCoords", "[Particle]")
{
  {
    constexpr auto mct = CoordsType::POS;
    auto mc_coords     = MCCoords<mct>();
    REQUIRE(mc_coords.positions.size() == 0);
    mc_coords.resize(3);
    REQUIRE(mc_coords.positions.size() == 3);
  }
  {
    constexpr auto mct = CoordsType::POS_SPIN;
    auto mc_coords     = MCCoords<mct>();
    REQUIRE(mc_coords.spins.size() == 0);
    mc_coords.resize(3);
    REQUIRE(mc_coords.positions.size() == 3);
    REQUIRE(mc_coords.spins.size() == 3);
  }
}

TEST_CASE("Taus", "[Particle]")
{
  MCCoords<CoordsType::POS> mc_coords_rs;
  auto tau     = 1.0;
  auto invmass = 0.2;
  auto taus_rs{makeTaus(mc_coords_rs, tau, invmass)};
  CHECK(Approx(taus_rs.tauovermass) == 0.2);
  CHECK(Approx(taus_rs.oneover2tau) == 2.5);
  CHECK(Approx(taus_rs.sqrttau) == 0.447213595499957927703605);
  MCCoords<CoordsType::POS_SPIN> mc_coords_rsspins;
  auto spin_mass    = 0.5;
  auto taus_rsspins = makeTaus(mc_coords_rsspins, tau, invmass, spin_mass);
  CHECK(Approx(taus_rsspins.spin_tauovermass) == 0.4);
  CHECK(Approx(taus_rsspins.spin_oneover2tau) == 1.25);
  CHECK(Approx(taus_rsspins.spin_sqrttau) == 0.632455532033675882352952);
}

} // namespace qmcplusplus
