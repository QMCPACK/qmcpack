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

} // namespace qmcplusplus
