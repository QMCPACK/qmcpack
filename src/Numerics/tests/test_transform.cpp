//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include <stdio.h>
#include <string>
#include <iostream>

#include "Numerics/Transform2GridFunctor.h"
#include "Numerics/OneDimQuinticSpline.h"

using std::string;

namespace qmcplusplus
{
class Input
{
public:
  typedef double real_type;

  double f(double r) { return r * r; }

  double df(double r) { return 2 * r; }
};

TEST_CASE("transform2gridfunctor", "[numerics]")
{
  typedef OneDimGridBase<double> GridType;
  typedef OneDimQuinticSpline<double> OutputType;

  GridType* agrid = new LogGrid<double>;
  agrid->set(0.1, 10, 10);
  OutputType output(agrid);
  Input input;
  Transform2GridFunctor<Input, OutputType> transform(input, output);
  double rmin = 0.1;
  double rmax = 10;
  int npts    = 10;
  transform.generate(rmin, rmax, npts);
  REQUIRE(output.splint(0.1) == Approx(0.01));
  REQUIRE(output.splint(0.15) == Approx(0.0225));
  REQUIRE(output.splint(7.0) == Approx(49.0));
  REQUIRE(output.splint(10) == Approx(100.0));
  delete agrid;
}
} // namespace qmcplusplus
