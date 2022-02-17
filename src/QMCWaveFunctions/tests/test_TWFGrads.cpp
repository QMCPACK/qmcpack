//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "TWFGrads.hpp"

namespace qmcplusplus
{

TEST_CASE("TWFGrads", "[QMCWaveFunctions]")
{
  {
    constexpr auto CT = CoordsType::POS;
    auto grads        = TWFGrads<CT>();
    REQUIRE(grads.grads_positions.size() == 0);
    grads.resize(7);
    REQUIRE(grads.grads_positions.size() == 7);
  }
  {
    constexpr auto CT = CoordsType::POS_SPIN;
    auto grads    = TWFGrads<CT>();
    REQUIRE(grads.grads_positions.size() == 0);
    REQUIRE(grads.grads_spins.size() == 0);
    grads.resize(5);
    REQUIRE(grads.grads_positions.size() == 5);
    REQUIRE(grads.grads_spins.size() == 5);
  }
}

} // namespace qmcplusplus
