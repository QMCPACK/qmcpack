//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"
#include "Numerics/SplineBound.hpp"

namespace qmcplusplus
{
template<typename T>
void test_spline_bounds()
{
  T x = 2.2;
  T dx;
  int ind;
  int ng = 10;
  getSplineBound(x, ng, ind, dx);
  CHECK(dx == Approx(0.2));
  REQUIRE(ind == 2);

  // check clamping to a maximum index value
  x = 10.5;
  getSplineBound(x, ng, ind, dx);
  CHECK(dx == Approx(0.5));
  REQUIRE(ind == 10);

  x = 11.5;
  getSplineBound(x, ng, ind, dx);
  CHECK(dx == Approx(1.0));
  REQUIRE(ind == 10);

  // check clamping to a zero index value
  x = -1.3;
  getSplineBound(x, ng, ind, dx);
  CHECK(dx == Approx(0.0));
  REQUIRE(ind == 0);
}

TEST_CASE("getSplineBound double", "[numerics]") { test_spline_bounds<double>(); }
TEST_CASE("getSplineBound float", "[numerics]") { test_spline_bounds<float>(); }
} // namespace qmcplusplus
