//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//                    Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
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
    auto mc_coords     = MCCoords<mct>(3);
    REQUIRE(mc_coords.positions.size() == 3);

    QMCTraits::PosType p1({0.1, 0.2, 0.3});
    QMCTraits::PosType p2({0.4, 0.5, 0.6});
    QMCTraits::PosType p3({0.7, 0.8, 0.9});
    QMCTraits::PosType p4({-1.0, -1.0, -1.0});

    mc_coords.positions = {p1, p2, p3};
    auto shift_coords   = MCCoords<mct>(3);
    shift_coords.positions = {p4, p4, p4};
    mc_coords += shift_coords;
    CHECK(Approx(mc_coords.positions[0][0]) == -0.9);
    CHECK(Approx(mc_coords.positions[0][1]) == -0.8);
    CHECK(Approx(mc_coords.positions[0][2]) == -0.7);
    CHECK(Approx(mc_coords.positions[1][0]) == -0.6);
    CHECK(Approx(mc_coords.positions[1][1]) == -0.5);
    CHECK(Approx(mc_coords.positions[1][2]) == -0.4);
    CHECK(Approx(mc_coords.positions[2][0]) == -0.3);
    CHECK(Approx(mc_coords.positions[2][1]) == -0.2);
    CHECK(Approx(mc_coords.positions[2][2]) == -0.1);

    auto shift_coords2   = MCCoords<mct>(3);
    shift_coords2.positions = {-p4, -p4, -p4};
    mc_coords += shift_coords2;
    CHECK(Approx(mc_coords.positions[0][0]) == 0.1);
    CHECK(Approx(mc_coords.positions[0][1]) == 0.2);
    CHECK(Approx(mc_coords.positions[0][2]) == 0.3);
    CHECK(Approx(mc_coords.positions[1][0]) == 0.4);
    CHECK(Approx(mc_coords.positions[1][1]) == 0.5);
    CHECK(Approx(mc_coords.positions[1][2]) == 0.6);
    CHECK(Approx(mc_coords.positions[2][0]) == 0.7);
    CHECK(Approx(mc_coords.positions[2][1]) == 0.8);
    CHECK(Approx(mc_coords.positions[2][2]) == 0.9);

    MCCoords<mct> subset(1);
    mc_coords.getSubset(1, 1, subset);
    REQUIRE(subset.positions.size() == 1);
    CHECK(Approx(subset.positions[0][0]) == 0.4);
    CHECK(Approx(subset.positions[0][1]) == 0.5);
    CHECK(Approx(subset.positions[0][2]) == 0.6);
  }
  {
    constexpr auto mct = CoordsType::POS_SPIN;
    auto mc_coords     = MCCoords<mct>(3);
    REQUIRE(mc_coords.positions.size() == 3);
    REQUIRE(mc_coords.spins.size() == 3);

    mc_coords.spins   = {0.1, 0.2, 0.3};
    auto shift_coords = MCCoords<mct>(3);
    shift_coords.spins = {1.0, 1.0, 1.0};
    mc_coords += shift_coords;
    CHECK(Approx(mc_coords.spins[0]) == 1.1);
    CHECK(Approx(mc_coords.spins[1]) == 1.2);
    CHECK(Approx(mc_coords.spins[2]) == 1.3);

    auto shift_coords2   = MCCoords<mct>(3);
    shift_coords2.spins = {-1.0, -1.0, -1.0};
    mc_coords += shift_coords2;
    CHECK(Approx(mc_coords.spins[0]) == 0.1);
    CHECK(Approx(mc_coords.spins[1]) == 0.2);
    CHECK(Approx(mc_coords.spins[2]) == 0.3);

    MCCoords<mct> subset(1);
    mc_coords.getSubset(2, 1, subset);
    REQUIRE(subset.spins.size() == 1);
    CHECK(Approx(subset.spins[0]) == 0.3);
  }
}

} // namespace qmcplusplus
