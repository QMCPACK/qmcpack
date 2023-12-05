//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "Numerics/OneDimCubicSpline.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
// Generated from gen_cubic_spline.py
TEST_CASE("spline_function_1", "[numerics]")
{
  const int n = 3;
  double x[n] = {0.0, 1.0, 2.0};
  double y[n] = {1.0, 2.0, 1.5};
  double y2[n];

  CubicSplineSolve(x, y, n, 1e+33, 1e+33, y2);

  CHECK(y2[0] == Approx(0));
  CHECK(y2[1] == Approx(-2.25));
  CHECK(y2[2] == Approx(0));
}


// Generated from gen_cubic_spline.py
TEST_CASE("spline_function_2", "[numerics]")
{
  const int n = 3;
  double x[n] = {0.0, 1.0, 2.0};
  double y[n] = {1.0, 2.0, 1.5};
  double y2[n];

  CubicSplineSolve(x, y, n, 1e+33, 1.0, y2);

  CHECK(y2[0] == Approx(0));
  CHECK(y2[1] == Approx(-3.85714));
  CHECK(y2[2] == Approx(6.42857));
}


// Generated from gen_cubic_spline.py
TEST_CASE("spline_function_3", "[numerics]")
{
  const int n = 3;
  double x[n] = {0.0, 1.0, 2.0};
  double y[n] = {1.0, 2.0, 1.5};
  double y2[n];

  CubicSplineSolve(x, y, n, 1.0, 2.0, y2);

  CHECK(y2[0] == Approx(2.75));
  CHECK(y2[1] == Approx(-5.5));
  CHECK(y2[2] == Approx(10.25));
}


// Generated from gen_cubic_spline.py
TEST_CASE("spline_function_4", "[numerics]")
{
  const int n = 4;
  double x[n] = {0.0, 1.2, 2.4, 3.0};
  double y[n] = {1.0, 2.0, 1.5, 1.8};
  double y2[n];

  CubicSplineSolve(x, y, n, 1e+33, 1e+33, y2);

  CHECK(y2[0] == Approx(0));
  CHECK(y2[1] == Approx(-2.12121));
  CHECK(y2[2] == Approx(2.23485));
  CHECK(y2[3] == Approx(0));
}


// Generated from gen_cubic_spline.py
TEST_CASE("spline_function_5", "[numerics]")
{
  const int n = 4;
  double x[n] = {0.0, 1.2, 2.4, 3.0};
  double y[n] = {1.0, 2.0, 1.5, 1.8};
  double y2[n];

  CubicSplineSolve(x, y, n, 0.0, 1.0, y2);

  CHECK(y2[0] == Approx(3.60507));
  CHECK(y2[1] == Approx(-3.04348));
  CHECK(y2[2] == Approx(2.31884));
  CHECK(y2[3] == Approx(1.34058));
}


// Structure for holding value, first and second derivatives
struct D2U
{
  D2U(double val_, double du_, double d2u_) : val(val_), du(du_), d2u(d2u_) {}
  double val;
  double du;
  double d2u;
};


// Generated from gen_cubic_spline.py
TEST_CASE("one_dim_cubic_spline_1", "[numerics]")
{
  const int n               = 3;
  std::vector<double> yvals = {1.0, 2.0, 1.5};

  auto grid = std::make_unique<LinearGrid<double>>();
  grid->set(0.0, 2.0, n);

  OneDimCubicSpline<double> cubic_spline(std::move(grid), yvals);

  int imin = 0;
  int imax = 3 - 1;

  double yp0 = 1.0;
  double ypn = 2.0;

  cubic_spline.spline(imin, yp0, imax, ypn);

  std::vector<double> check_xvals = {0.0, 0.39999999998, 0.79999999996, 1.19999999994, 1.59999999992, 1.9999999999};
  std::vector<double> check_yvals;
  std::vector<D2U> check_yvals_d2u;

  for (int i = 0; i < check_xvals.size(); i++)
  {
    double r   = check_xvals[i];
    double val = cubic_spline.splint(r);
    check_yvals.push_back(val);

    double du, d2u;
    double val2 = cubic_spline.splint(r, du, d2u);
    check_yvals_d2u.push_back(D2U(val2, du, d2u));

    //std::cout << i << " r = " << r << " val = " << val << " " << check_yvals[i] << std::endl;
  }

  CHECK(check_yvals[0] == Approx(1));
  CHECK(check_yvals[1] == Approx(1.532));
  CHECK(check_yvals[2] == Approx(1.976));
  CHECK(check_yvals[3] == Approx(1.836));
  CHECK(check_yvals[4] == Approx(1.352));
  CHECK(check_yvals[5] == Approx(1.5));

  CHECK(check_yvals_d2u[0].val == Approx(1));
  CHECK(check_yvals_d2u[0].du == Approx(1));
  CHECK(check_yvals_d2u[0].d2u == Approx(2.75));
  CHECK(check_yvals_d2u[1].val == Approx(1.532));
  CHECK(check_yvals_d2u[1].du == Approx(1.44));
  CHECK(check_yvals_d2u[1].d2u == Approx(-0.55));
  CHECK(check_yvals_d2u[2].val == Approx(1.976));
  CHECK(check_yvals_d2u[2].du == Approx(0.56));
  CHECK(check_yvals_d2u[2].d2u == Approx(-3.85));
  CHECK(check_yvals_d2u[3].val == Approx(1.836));
  CHECK(check_yvals_d2u[3].du == Approx(-1.16));
  CHECK(check_yvals_d2u[3].d2u == Approx(-2.35));
  CHECK(check_yvals_d2u[4].val == Approx(1.352));
  CHECK(check_yvals_d2u[4].du == Approx(-0.84));
  CHECK(check_yvals_d2u[4].d2u == Approx(3.95));
  CHECK(check_yvals_d2u[5].val == Approx(1.5));
  CHECK(check_yvals_d2u[5].du == Approx(2));
  CHECK(check_yvals_d2u[5].d2u == Approx(10.25));
}

} // namespace qmcplusplus
