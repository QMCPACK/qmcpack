//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include <cmath>
#include "CPU/math.hpp"

namespace qmcplusplus
{

template<typename T>
void test_isnan()
{
  const T one(1);
  T a = std::sqrt(-one);
  T b = -a;
  CHECK(!qmcplusplus::isnan(one));
  CHECK(qmcplusplus::isnan(a));
  CHECK(qmcplusplus::isnan(b));
}

TEST_CASE("isnan", "[numerics]")
{
  test_isnan<float>();
  test_isnan<double>();
}

} // namespace qmcplusplus
