//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Numerics/Bessel.h"

namespace qmcplusplus
{

TEST_CASE("Bessel", "[numerics]")
{
  double bessel_array[100];
  bessel_steed_array_cpu(10, 1.0, bessel_array);

  REQUIRE(bessel_array[0] == Approx(0.84147098480789650670));
  REQUIRE(bessel_array[1] == Approx(0.30116867893975678925));
  REQUIRE(bessel_array[3] == Approx(9.006581117112515480e-03));
  REQUIRE(bessel_array[5] == Approx(9.256115861125815763e-05));
  REQUIRE(bessel_array[7] == Approx(4.790134198739488623e-07));
  REQUIRE(bessel_array[10] == Approx(7.116552640047313024e-11));

}

}
