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
#include "SizeLimitedDataQueue.hpp"

namespace qmcplusplus
{

TEST_CASE("SizeLimitedDataQueue", "[estimators]")
{
  SizeLimitedDataQueue<double, 1> weight_and_energy(3);
  CHECK(weight_and_energy.size() == 0);
  {
    weight_and_energy.push({1.0, {2.0}});
    CHECK(weight_and_energy.size() == 1);
    auto avg = weight_and_energy.weighted_avg();
    CHECK(Approx(avg[0]) == 2.0);
  }
  {
    weight_and_energy.push({3.0, {1.0}});
    CHECK(weight_and_energy.size() == 2);
    auto avg = weight_and_energy.weighted_avg();
    CHECK(Approx(avg[0]) == 1.25);
  }
  {
    SizeLimitedDataQueue<double, 1>::HistoryElement temp{0.5, {3.0}};
    weight_and_energy.push(std::move(temp));
    CHECK(weight_and_energy.size() == 3);
    auto avg = weight_and_energy.weighted_avg();
    CHECK(Approx(avg[0]) == 1.444444444);
  }
  {
    weight_and_energy.push({0.5, {3.0}});
    CHECK(weight_and_energy.size() == 3);
    auto avg = weight_and_energy.weighted_avg();
    CHECK(Approx(avg[0]) == 1.5);
  }
}

} // namespace qmcplusplus
