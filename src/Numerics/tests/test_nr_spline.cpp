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


#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "Numerics/NRSplineFunctions.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{

TEST_CASE("NR_spline_functions", "[numerics]")
{

  const int n = 2;
  double x[n] = {0.0, 1.0};
  double y[n] = {1.0, 2.0};
  double y2[n];

  NRCubicSpline(x, y, n, 1.0, 2.0, y2);

  REQUIRE(y2[0] != 0.0);
}

}

