//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Utilities/RandomGenerator.h"
#include "Utilities/BoostRandom.h"
#include <stdio.h>
#include <string>
#include <vector>

using std::string;
using std::vector;
namespace qmcplusplus
{
#ifdef HAVE_LIBBOOST
TEST_CASE("boost", "[utilities]")
{
  double d = Random();
  REQUIRE(d >= 0.0);
  REQUIRE(d < 1.0);
}

TEST_CASE("boost_mt19937_determinism", "[utilities]")
{
  // Verify BoostRandom (MT19937 internally) generates a fixed sequence given an fixed initial seed
  // This is not guaranteed by Boost but appears to be the case
  BoostRandom<OHMMS_PRECISION_FULL> our_rng(13);
  std::vector<OHMMS_PRECISION_FULL> expected = {0.7777024102, 0.6073413305, 0.237541216};
  for (auto i = 0; i < expected.size(); ++i)
  {
    REQUIRE(our_rng() == Approx(expected[i]));
  }
}

#endif

TEST_CASE("make_seed", "[utilities]")
{
  // not sure what to test here - mostly that it doesn't crash
  // It's based on 'time' so it will return different values at different calls
  // If the time is the same (1 second resolution), different inputs should
  //  give different seeds.
  unsigned int seed1 = make_seed(0, 0);
  unsigned int seed2 = make_seed(1, 1);

  REQUIRE(seed1 != seed2);
}

} // namespace qmcplusplus
