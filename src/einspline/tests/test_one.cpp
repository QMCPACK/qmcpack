//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "einspline/bspline_base.h"
#include "einspline/bspline_create.h"
#include "einspline/bspline_eval_d.h"
#include "einspline/multi_bspline_create.h"
#include "einspline/multi_bspline_eval_d.h"
#include "einspline/multi_bspline_eval_s.h"

#include <stdio.h>
#include <vector>
#include "config/stdlib/Constants.h"

TEST_CASE("double_1d_natural", "[einspline]")
{
  // two point case
  Ugrid x_grid;
  x_grid.start = 1.0;
  x_grid.end = 10.0;
  x_grid.num = 2;

  std::vector<double> data(2);
  data[0] = 2.0; data[1] = 3.0;

  BCtype_d xBC;
  xBC.lCode = NATURAL;
  xBC.rCode = NATURAL;
  auto s =
      std::unique_ptr<UBspline_1d_d, void (*)(void*)>{create_UBspline_1d_d(x_grid, xBC, data.data()), destroy_Bspline};

  REQUIRE(s);

  double val;

  eval_UBspline_1d_d(s.get(), 1.0, &val);
  CHECK(val == Approx(2.0));

  eval_UBspline_1d_d(s.get(), 9.9999999, &val);
  CHECK(val == Approx(3.0));

  // This should assert
  // eval_UBspline_1d_d(s.get(), 10.0, &val);
  // CHECK(val == Approx(3.0));

  eval_UBspline_1d_d(s.get(), 5.5, &val);
  CHECK(val == Approx(2.5));


  // three point case
  x_grid.start = 1.0;
  x_grid.end = 10.0;
  x_grid.num = 3;

  data.resize(3);
  data[0] = 2.0; data[1] = 2.7; data[2] = 3.0;

  xBC.lCode = NATURAL;
  xBC.rCode = NATURAL;
  s.reset(create_UBspline_1d_d(x_grid, xBC, data.data()));

  REQUIRE(s);

  eval_UBspline_1d_d(s.get(), 1.0, &val);
  CHECK(val == Approx(2.0));

  eval_UBspline_1d_d(s.get(), 9.9999999, &val);
  CHECK(val == Approx(3.0));

  eval_UBspline_1d_d(s.get(), 5.5, &val);
  CHECK(val == Approx(2.7));
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
  auto s    = std::unique_ptr<multi_UBspline_1d_d, void (*)(void*)>{create_multi_UBspline_1d_d(x_grid, xBC, 1),
                                                                 destroy_Bspline};
  REQUIRE(s);

  set_multi_UBspline_1d_d(s.get(), 0, data);

  double val;
  eval_multi_UBspline_1d_d(s.get(), 1.0, &val);
  CHECK(val == Approx(2.0));
}

TEST_CASE("double_1d_periodic", "[einspline]")
{
  Ugrid x_grid;
  x_grid.start = 0.0;
  x_grid.end = 1.0;
  //Enough grid points are required to do the micro evaluation test.
  constexpr int N = 12;
  x_grid.num = N;
  double delta = (x_grid.end - x_grid.start)/x_grid.num;

  double tpi = 2*M_PI;
  double data[N];
  for (int i = 0; i < N; i++)
  {
    double x = delta*i;
    data[i] = sin(tpi*x);
  }

  BCtype_d bc;
  bc.lCode = PERIODIC;
  bc.rCode = PERIODIC;

  auto s = std::unique_ptr<UBspline_1d_d, void (*)(void*)>{create_UBspline_1d_d(x_grid, bc, data), destroy_Bspline};

  REQUIRE(s);

  double val;
  eval_UBspline_1d_d(s.get(), 0.0, &val);
  CHECK(val == Approx(0.0));

  eval_UBspline_1d_d(s.get(), delta, &val);
  CHECK(val == Approx(data[1]));

  double micro_delta = delta / 4.0;
  int micro_N = N * 4;
  double micro_data[N*4];
  for (int i = 0 ; i < micro_N; i++)
  {
    double x = micro_delta * i;
    micro_data[i] = sin(tpi*x);
  }
  eval_UBspline_1d_d(s.get(), micro_delta * 3, &val);
  CHECK(val == Approx(micro_data[3]).epsilon(0.001));

  eval_UBspline_1d_d(s.get(), micro_delta * 17, &val);
  CHECK(val == Approx(micro_data[17]).epsilon(0.001));

  eval_UBspline_1d_d(s.get(), micro_delta * 31, &val);
  CHECK(val == Approx(micro_data[31]).epsilon(0.001));
}
