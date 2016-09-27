#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana Champain
#//
#// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana Champain
#//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "einspline/bspline_base.h"
#include "einspline/bspline_create.h"
#include "einspline/bspline_eval_d.h"

#include <stdio.h>
#include <string>

using std::string;

TEST_CASE("double_3d_natural","[einspline]")
{
  Ugrid grid;
  grid.start = 1.0;
  grid.end = 10.0;
  grid.num = 2;

  double data[8];
  data[0] = 2.0;
  data[1] = 3.0;
  data[2] = 4.0;
  data[3] = 5.0;
  data[4] = 6.0;
  data[5] = 7.0;
  data[6] = 8.0;
  data[7] = 9.0;

  BCtype_d bc;
  bc.lCode = NATURAL;
  bc.rCode = NATURAL;
  UBspline_3d_d* s = create_UBspline_3d_d(grid, grid, grid, bc, bc, bc, data);

  REQUIRE(s);

  double val;
  eval_UBspline_3d_d(s, 1.0, 1.0, 1.0, &val);
  REQUIRE(val == Approx(2.0));

  eval_UBspline_3d_d(s, 10.0, 10.0, 10.0, &val);
  REQUIRE(val == Approx(9.0));

  destroy_Bspline(s);
}
