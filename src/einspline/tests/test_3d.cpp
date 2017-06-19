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


#include "catch.hpp"

#include "einspline/bspline_base.h"
#include "einspline/bspline_create.h"
#include "einspline/bspline_eval_d.h"
#include "einspline/multi_bspline_eval_d.h"
#include "einspline/multi_bspline_eval_s.h"
#include "einspline/multi_bspline_create.h"

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

TEST_CASE("double_3d_periodic","[einspline]")
{
  Ugrid grid;
  grid.start = 0.0;
  grid.end = 1.0;
  int N = 5;
  grid.num = N;
  double delta = (grid.end - grid.start)/grid.num;

  double tpi = 2*M_PI;
  double data[N*N*N];
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      for (int k = 0; k < N; k++) {
        double x = delta*i;
        double y = delta*j;
        double z = delta*k;
        data[N*N*i + N*j + k] = sin(tpi*x) + sin(3*tpi*y) + sin(4*tpi*z);
      }
    }
  }

  BCtype_d bc;
  bc.lCode = PERIODIC;
  bc.rCode = PERIODIC;

  UBspline_3d_d* s = create_UBspline_3d_d(grid, grid, grid, bc, bc, bc, data);
  REQUIRE(s);

  double val;
  eval_UBspline_3d_d(s, 0.0, 0.0, 0.0, &val);
  REQUIRE(val == Approx(0.0));

  double x = delta;
  eval_UBspline_3d_d(s, delta, delta, delta, &val);
  REQUIRE(val == Approx(data[N*N + N + 1]));

  multi_UBspline_3d_d* m = create_multi_UBspline_3d_d(grid, grid, grid, bc, bc, bc, 1);
  REQUIRE(m);

  set_multi_UBspline_3d_d(m, 0, data);

  eval_multi_UBspline_3d_d(m, delta, delta, delta, &val);
  REQUIRE(val == Approx(data[N*N + N + 1]));


  BCtype_s bc_s;
  bc_s.lCode = PERIODIC;
  bc_s.rCode = PERIODIC;

  multi_UBspline_3d_s* ms = create_multi_UBspline_3d_s(grid, grid, grid, bc_s, bc_s, bc_s, 1);
  REQUIRE(ms);
  set_multi_UBspline_3d_s_d(ms, 0, data);

  float fval;
  eval_multi_UBspline_3d_s(ms, delta, delta, delta, &fval);
  REQUIRE(fval == Approx(data[N*N + N + 1]));

  float grads[3];
  float hess[9];
  eval_multi_UBspline_3d_s_vgh(ms, 0.1, 0.2, 0.0, &fval, grads, hess);

  // See miniqmc-python/splines/test_3d.py for values
  REQUIRE(grads[0] == Approx(5.11104213833));
  REQUIRE(grads[1] == Approx(5.98910634201));
  REQUIRE(grads[2] == Approx(-6.17832080889));

  // All off-diagonal values of the Hessian for this data should be zero
  REQUIRE(hess[1] == Approx(0.0));
  REQUIRE(hess[2] == Approx(0.0));

  REQUIRE(hess[3] == Approx(0.0));
  REQUIRE(hess[5] == Approx(0.0));

  REQUIRE(hess[6] == Approx(0.0));
  REQUIRE(hess[7] == Approx(0.0));
}
