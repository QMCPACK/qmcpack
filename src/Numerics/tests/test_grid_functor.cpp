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

#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"


#include <stdio.h>
#include <string>
#include <sstream>

using std::string;
using std::stringstream;

namespace qmcplusplus
{

TEST_CASE("double_1d_grid_functor", "[numerics]")
{
  LinearGrid<double> grid;
  OneDimGridFunctor<double> f(&grid);

  grid.set(0.0, 1.0, 3);

  REQUIRE(grid.size() == 3);
  REQUIRE(grid.rmin() == 0.0);
  REQUIRE(grid.rmax() == 1.0);
  REQUIRE(grid.dh() == Approx(0.5));
  REQUIRE(grid.dr(1) == Approx(0.5));

}

TEST_CASE("double_1d_grid_functor_vs_n", "[numerics]")
{
  for (int n = 2; n < 5; n++) {
    stringstream sec_name;
    sec_name << "grid size " << n;
    SECTION(sec_name.str())
    {
      LinearGrid<double> grid;
      grid.set(0.0, 1.0, n);
      REQUIRE(grid.size() == n);
      REQUIRE(grid.rmin() == 0.0);
      REQUIRE(grid.rmax() == 1.0);
      REQUIRE(grid.dh() == Approx(1.0/(n-1)));
      REQUIRE(grid.dr(0) == Approx(1.0/(n-1)));
    }
  }
}

}
