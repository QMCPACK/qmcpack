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
  // two point case
  Ugrid x_grid;
  x_grid.start = 1.0;
  x_grid.end = 10.0;
  x_grid.num = 2;

  double data[3];
  data[0] = 2.0; data[1] = 3.0;

  BCtype_d xBC;
  xBC.lCode = NATURAL;
  xBC.rCode = NATURAL;
  UBspline_1d_d* s = create_UBspline_1d_d(x_grid, xBC, data);

  REQUIRE(s);

  double val;
  eval_UBspline_1d_d(s, 1.0, &val);
  REQUIRE(val == Approx(2.0));

  eval_UBspline_1d_d(s, 9.9999999, &val);
  REQUIRE(val == Approx(3.0));

  // NO! this results in ipart = 1 reads outside of the coefs allocation
  // No code that ever calls this function with xBC.lCode = NATURAL may have x >= xgrid.end
  eval_UBspline_1d_d(s, 10.0, &val);
  REQUIRE(val == Approx(3.0));

  eval_UBspline_1d_d(s, 5.5, &val);
  REQUIRE(val == Approx(2.5));

  destroy_Bspline(s);

  // three point case
  x_grid.start = 1.0;
  x_grid.end = 10.0;
  x_grid.num = 3;

  data[0] = 2.0; data[1] = 2.7; data[2] = 3.0;

  xBC.lCode = NATURAL;
  xBC.rCode = NATURAL;
  s = create_UBspline_1d_d(x_grid, xBC, data);

  REQUIRE(s);

  eval_UBspline_1d_d(s, 1.0, &val);
  REQUIRE(val == Approx(2.0));

  eval_UBspline_1d_d(s, 9.9999999, &val);
  REQUIRE(val == Approx(3.0));

  // NO! this results in ipart = 1 reads outside of the coefs allocation
  // No code that ever calls this function with xBC.lCode = NATURAL may have x >= xgrid.end
  eval_UBspline_1d_d(s, 10.0, &val);
  REQUIRE(val == Approx(3.0));

  eval_UBspline_1d_d(s, 5.5, &val);
  REQUIRE(val == Approx(2.7));

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
void test_multi(multi_UBspline_1d_s *cpuSpline, float *pos, float *vals_cuda);

TEST_CASE("multi_cuda_wrapper", "[einspline]")
{
  Ugrid x_grid;
  // GPU versions require the grid to start at zero
  x_grid.start = 0.0;
  x_grid.end = 10.0;
  x_grid.num = 2;

  float data[2];
  data[0] = 2.0;
  data[1] = 3.0;

  BCtype_s xBC;
  xBC.lCode = NATURAL;
  xBC.rCode = NATURAL;
  multi_UBspline_1d_s* s = create_multi_UBspline_1d_s(x_grid, xBC, 1);
  REQUIRE(s);
  set_multi_UBspline_1d_s(s, 0, data);

  float pos[1];
  pos[0] = 0.0;

  // Check the CPU value
  float cpu_val[1];
  eval_multi_UBspline_1d_s(s, pos[0], cpu_val);
  REQUIRE(cpu_val[0] == 2.0);

  pos[0] = 11.0;
  // Check the CPU value
  eval_multi_UBspline_1d_s(s, pos[0], cpu_val);
  REQUIRE(cpu_val[0] == 3.0);

  // Check the GPU value
  float vals_output[1];
  test_multi(s, pos, vals_output);
  REQUIRE(vals_output[0] == 2.0);
}

#endif
