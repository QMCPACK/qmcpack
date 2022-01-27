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


TEST_CASE("MCCoords", "[QMCDrivers]")
{
  {
    constexpr auto mct = CoordsTypes::RS;
    auto mc_coords     = MCCoords<mct>();
    REQUIRE(mc_coords.positions.size() == 0);
    mc_coords.resize(3);
    REQUIRE(mc_coords.positions.size() == 3);
  }
  {
    constexpr auto mct = CoordsTypes::RSSPINS;
    auto mc_coords     = MCCoords<mct>();
    REQUIRE(mc_coords.spins.size() == 0);
    mc_coords.resize(3);
    REQUIRE(mc_coords.positions.size() == 3);
    REQUIRE(mc_coords.spins.size() == 3);
  }
}

TEST_CASE("Taus", "[QMCDrivers]")
{
  MCCoords<CoordsTypes::RS> mc_coords_rs;
  auto tau     = 1.0;
  auto invmass = 0.2;
  auto taus_rs{makeTaus(mc_coords_rs, tau, invmass)};
  CHECK(taus_rs.tauovermass == 0.2);
  CHECK(taus_rs.oneover2tau == 2.5);
  CHECK(taus_rs.sqrttau == 0.447213595499957927703605);
  MCCoords<CoordsTypes::RSSPINS> mc_coords_rsspins;
  auto spin_mass    = 0.5;
  auto taus_rsspins = makeTaus(mc_coords_rsspins, tau, invmass, spin_mass);
  CHECK(taus_rsspins.spin_tauovermass == 0.4);
  CHECK(taus_rsspins.spin_oneover2tau == 1.25);
  CHECK(taus_rsspins.spin_sqrttau == 0.632455532033675882352952);
}

TEST_CASE("generateDeltas", "[QMCDrivers]")
{
  int size_test = 7;
  std::vector<double> gauss_random_vals(size_test * 3 + (size_test * 3) % 2 + size_test );
  {
    StdRandom<double> rng;
    makeGaussRandomWithEngine(gauss_random_vals, rng);
  }

  auto checkRs = [&](auto& rs) {
    for (int i = 0; i < size_test; ++i)
      CHECK(TinyVector<double, 3>(gauss_random_vals[3 * i], gauss_random_vals[3 * i + 1],
                                  gauss_random_vals[3 * i + 2]) == rs[i]);
  };

  MCCoords<CoordsTypes::RS> mc_coords_rs;
  mc_coords_rs.resize(size_test);
  {
    StdRandom<double> rng;
    makeGaussRandomWithEngine(mc_coords_rs, rng);
    std::cout << mc_coords_rs.positions;
    checkRs(mc_coords_rs.positions);
  }
  MCCoords<CoordsTypes::RSSPINS> mc_coords_rsspins;
  mc_coords_rsspins.resize(size_test);
  {
    StdRandom<double> rng;
    makeGaussRandomWithEngine(mc_coords_rsspins, rng);
    checkRs(mc_coords_rsspins.positions);
    // Mod 2 is result of how gaussianDistribution is generated.
    int offset_for_rs = ( 3 * size_test ) + (3* size_test) % 2;
    for (int i = 0; i < size_test; ++i)
      CHECK(typename decltype(mc_coords_rsspins.spins)::value_type{gauss_random_vals[offset_for_rs + i]} ==
            mc_coords_rsspins.spins[i]);
  }
}

} // namespace qmcplusplus
