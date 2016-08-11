
#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "Numerics/NRSplineFunctions.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{

TEST_CASE("NR_spline_functions", "[numerics]")
{

  const int n = 2;
  double x[n] = {0.0, 1.0};
  double y[n] = {1.0, 2.0};
  double y2[n];

  NRCubicSpline(x, y, n, 1.0, 2.0, y2);

  REQUIRE(y2[0] != 0.0);
}

}

