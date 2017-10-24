#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//
#// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////


#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "einspline/bspline_base.h"
#include "einspline/bspline_create.h"
#include "einspline/bspline_eval_d.h"
#include "einspline/multi_bspline_create.h"
#include "einspline/multi_bspline_eval_d.h"
#include "einspline/multi_bspline_eval_s.h"

#include <stdio.h>
#include <string>

using std::string;

TEST_CASE("double_1d_natural", "[einspline]")
{
  Ugrid x_grid;
  x_grid.start = 1.0;
  x_grid.end = 10.0;
  x_grid.num = 2;

  double data[2];
  data[0] = 2.0; data[1] = 3.0;

  BCtype_d xBC;
  xBC.lCode = NATURAL;
  xBC.rCode = NATURAL;
  UBspline_1d_d* s = create_UBspline_1d_d(x_grid, xBC, data);

  REQUIRE(s);

  double val;
  eval_UBspline_1d_d(s, 1.0, &val);
  REQUIRE(val == Approx(2.0));

  eval_UBspline_1d_d(s, 10.0, &val);
  REQUIRE(val == Approx(3.0));

  eval_UBspline_1d_d(s, 5.5, &val);
  REQUIRE(val == Approx(2.5));

  destroy_Bspline(s);
}

TEST_CASE("double_1d_multi", "[einspline]")
{
  Ugrid x_grid;
  x_grid.start = 1.0;
  x_grid.end = 10.0;
  x_grid.num = 2;

  double data[2];
  data[0] = 2.0;
  data[1] = 3.0;

  BCtype_d xBC;
  xBC.lCode = NATURAL;
  xBC.rCode = NATURAL;
  multi_UBspline_1d_d* s = create_multi_UBspline_1d_d(x_grid, xBC, 1);
  REQUIRE(s);

  set_multi_UBspline_1d_d(s, 0, data);

  double val;
  eval_multi_UBspline_1d_d(s, 1.0, &val);
  REQUIRE(val == Approx(2.0));
}

#ifdef QMC_CUDA
// If code from the .cu file is not called, the tests there don't get run.
// Call a simple function to force the link.
int force_cuda_link();

TEST_CASE("force_cuda_link", "[einspline]")
{
    int a = force_cuda_link();
}
#endif
