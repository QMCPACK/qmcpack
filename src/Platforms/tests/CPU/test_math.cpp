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
#include <catch2/catch_test_macros.hpp>

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

template<typename T>
void test_numeric_limits()
{
  const T inf = qmcplusplus::numeric_limits<T>::infinity();
  const T nan = qmcplusplus::numeric_limits<T>::quiet_NaN();

  CHECK(qmcplusplus::isinf(inf));
  CHECK(!qmcplusplus::isfinite(inf));
  CHECK(!qmcplusplus::isnan(inf));
  CHECK(inf > T(0));

  CHECK(qmcplusplus::isinf(-inf));
  CHECK(!qmcplusplus::isfinite(-inf));
  CHECK(!qmcplusplus::isnan(-inf));
  CHECK(-inf < T(0));

  CHECK(qmcplusplus::isnan(nan));
  CHECK(!qmcplusplus::isfinite(nan));
  CHECK(!qmcplusplus::isinf(nan));

  const T zero(0);
  CHECK(qmcplusplus::isfinite(zero));
  CHECK(!qmcplusplus::isnan(zero));
  CHECK(!qmcplusplus::isinf(zero));

  const T one(1);
  CHECK(qmcplusplus::isfinite(one));
  CHECK(!qmcplusplus::isnan(one));
  CHECK(!qmcplusplus::isinf(one));
}

TEST_CASE("numeric_limits_infinity_and_nan", "[numerics]")
{
  test_numeric_limits<float>();
  test_numeric_limits<double>();
}

} // namespace qmcplusplus
