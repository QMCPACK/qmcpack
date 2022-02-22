//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "Numerics/OneDimCubicSplineLinearGrid.h"

namespace qmcplusplus
{

// Copy and modified from one_dim_cubic_spline_1
TEST_CASE("test oneDimCubicSplineLinearGrid", "[numerics]")
{
  const int n               = 3;
  std::vector<double> yvals = {1.0, 2.0, 1.5};

  auto grid = std::make_unique<LinearGrid<double>>();
  grid->set(0.5, 2.0, n);

  OneDimCubicSpline<double> cubic_spline(std::move(grid), yvals);

  int imin = 0;
  int imax = 3 - 1;

  double yp0 = 1.0;
  double ypn = 2.0;

  cubic_spline.spline(imin, yp0, imax, ypn);

  OneDimCubicSplineLinearGrid linear_grid_cubic_spline(cubic_spline);

  std::vector<double> check_xvals = {0.0, 0.39999999998, 0.79999999996, 1.19999999994, 1.59999999992, 1.9999999999, 2.5};
  std::vector<double> check_yvals;

  for (int i = 0; i < check_xvals.size(); i++)
  {
    double r   = check_xvals[i];
    double val = cubic_spline.splint(r);
    check_yvals.push_back(val);

    //std::cout << i << " r = " << r << " val = " << val << " " << check_yvals[i] << std::endl;
  }

  CHECK(check_yvals[0] == Approx(0.5));
  CHECK(check_yvals[1] == Approx(0.9));
  CHECK(check_yvals[2] == Approx(1.4779999999));
  CHECK(check_yvals[3] == Approx(2.0012592592));
  CHECK(check_yvals[4] == Approx(1.575851852));
  CHECK(check_yvals[5] == Approx(1.5));
  CHECK(check_yvals[6] == Approx(1.5));
}

} // namespace qmcplusplus
