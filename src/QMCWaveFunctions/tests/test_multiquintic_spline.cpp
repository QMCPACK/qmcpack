
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

#include "Configuration.h"

#include "Message/Communicate.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Numerics/OneDimQuinticSpline.h"
#include "QMCWaveFunctions/LCAO/MultiQuinticSpline1D.h"

namespace qmcplusplus
{
TEST_CASE("LogGridLight", "[wavefunction][LCAO]")
{
  LogGridLight<double> grid;
  grid.set(0.1, 1.0, 5);

  int idx = grid.locate(0.1);
  REQUIRE(idx == 0);

  int idx2 = grid.locate(0.9);
  REQUIRE(idx2 == 3);

  double t = grid(0);
  CHECK(t == Approx(0.1));

  int idx3  = 0;
  double t2 = grid.getCLForQuintic(0.1, idx3);
  REQUIRE(idx3 == 0);
  CHECK(t2 == Approx(0.0));
}


TEST_CASE("MultiQuinticSpline", "[wavefunction][LCAO]")
{
  Vector<double> data(5);
  data[0] = 0.0;
  data[1] = 1.0;
  data[2] = 2.0;
  data[3] = 3.0;
  data[4] = 4.0;

  Vector<double> data2(5);
  data[0] = 2.0;
  data[1] = 9.0;
  data[2] = 1.0;
  data[3] = 2.0;
  data[4] = 0.1;

  auto agrid = std::make_unique<LogGrid<double>>();
  agrid->set(.1, 1.0, 5);

  OneDimQuinticSpline<double> spline1(agrid->makeClone());
  spline1.set(data);
  spline1.spline();

  OneDimQuinticSpline<double> spline2(agrid->makeClone());
  spline2.set(data2);
  spline2.spline();


  MultiQuinticSpline1D<double> m_spline;

  m_spline.initialize(*agrid, 2);
  m_spline.add_spline(0, spline1);
  m_spline.add_spline(1, spline2);


  // Compare the values from the multi-spline routines against the values from
  // the original single-spline routines.

  double u[2];
  double du[2];
  double d2u[2];
  double d3u[2];

  double u1, u2;
  double du1 = 0.0, du2 = 0.0;
  double d2u1 = 0.0, d2u2 = 0.0;
  double d3u1 = 0.0, d3u2 = 0.0;
  for (int i = 0; i < 10; i++)
  {
    double r = 0.08 * i + 0.1;
    m_spline.evaluate(r, u, du, d2u);

    u1 = spline1.splint(r, du1, d2u1);
    u2 = spline2.splint(r, du2, d2u2);

    CHECK(u[0] == Approx(u1));
    CHECK(du[0] == Approx(du1));
    CHECK(d2u[0] == Approx(d2u1));
    CHECK(u[1] == Approx(u2));
    CHECK(du[1] == Approx(du2));
    CHECK(d2u[1] == Approx(d2u2));

    m_spline.evaluate(r, u, du, d2u, d3u);
    u1 = spline1.splint(r, du1, d2u1, d3u1);
    u2 = spline2.splint(r, du2, d2u2, d3u2);

    CHECK(u[0] == Approx(u1));
    CHECK(du[0] == Approx(du1));
    CHECK(d2u[0] == Approx(d2u1));
    CHECK(d3u[0] == Approx(d3u1));
    CHECK(u[1] == Approx(u2));
    CHECK(du[1] == Approx(du2));
    CHECK(d2u[1] == Approx(d2u2));
    CHECK(d3u[1] == Approx(d3u2));
  }
}

} // namespace qmcplusplus
