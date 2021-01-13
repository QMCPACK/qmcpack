//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "RandomForTest.h"
#include <algorithm>

namespace qmcplusplus
{
namespace testing
{
RandomForTest::RandomForTest() { rng.init(0, 1, 111); }

std::vector<RandomForTest::Real> RandomForTest::getRealRandoms(int ncount)
{
  std::vector<Real> rng_reals;
  rng_reals.reserve(ncount);
  std::generate_n(std::back_inserter(rng_reals), ncount, rng);
  return rng_reals;
}

void RandomForTest::makeRngReals(std::vector<Real>& rngReals)
{
  // until c++ std = 17
  //std::generate(rng_reals.begin(), rng_reals.end(), rng());
  for (auto& rng_real : rngReals)
    rng_real = rng();
}

} // namespace testing
} // namespace qmcplusplus
