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
#include <catch2/catch_test_macros.hpp>
#include "Utilities/for_testing/Catch2Approx.h"

#include "Numerics/Bessel.h"

namespace qmcplusplus
{
TEST_CASE("Bessel", "[numerics]")
{
  double bessel_array[100];
  bessel_steed_array_cpu(10, 1.0, bessel_array);

  CHECK(bessel_array[0] == Approx(0.84147098480789650670));
  CHECK(bessel_array[1] == Approx(0.30116867893975678925));
  CHECK(bessel_array[3] == Approx(9.006581117112515480e-03));
  CHECK(bessel_array[5] == Approx(9.256115861125815763e-05));
  CHECK(bessel_array[7] == Approx(4.790134198739488623e-07));
  CHECK(bessel_array[10] == Approx(7.116552640047313024e-11));
}

TEST_CASE("Bessel edge and invalid inputs", "[numerics]")
{
  double zero_array[3];
  bessel_steed_array_cpu(2, 0.0, zero_array);
  CHECK(zero_array[0] == Approx(1.0));
  CHECK(zero_array[1] == Approx(0.0));
  CHECK(zero_array[2] == Approx(0.0));

  double tiny_array[3];
  bessel_steed_array_cpu(2, 1.0e-6, tiny_array);
  CHECK(tiny_array[0] == Approx(1.0).margin(1.0e-12));
  CHECK(tiny_array[1] == Approx(1.0e-6 / 3.0).margin(1.0e-12));
  CHECK(tiny_array[2] == Approx(0.0).margin(1.0e-18));

  double invalid_array[1];
  REQUIRE_THROWS_AS(bessel_steed_array_cpu(-1, 1.0, invalid_array), std::runtime_error);
  REQUIRE_THROWS_AS(bessel_steed_array_cpu(2, -1.0, invalid_array), std::runtime_error);
}

} // namespace qmcplusplus
