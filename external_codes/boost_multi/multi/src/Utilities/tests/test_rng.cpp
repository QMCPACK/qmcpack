//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by:  Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Utilities/RandomGenerator.h"
#include <stdio.h>
#include <string>
#include <vector>

using std::string;
using std::vector;
namespace qmcplusplus
{
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
