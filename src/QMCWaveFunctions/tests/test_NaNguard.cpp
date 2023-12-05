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
#include <NaNguard.h>
#include <cmath>

namespace qmcplusplus
{

TEST_CASE("NaNguard", "[wavefunction]")
{
  const QMCTraits::ValueType const_nan = std::sqrt(-1.0);
  CHECK_THROWS_WITH(NaNguard::checkOneParticleRatio(const_nan, "unit test"), Catch::Matchers::Contains("NaNguard::checkOneParticleRatio"));
  CHECK_NOTHROW(NaNguard::checkOneParticleRatio(0.0, "unit test"));
  CHECK_THROWS_WITH(NaNguard::checkOneParticleGradients({const_nan, -const_nan, const_nan}, "unit test"), Catch::Matchers::Contains("NaNguard::checkOneParticleGradients"));
  CHECK_NOTHROW(NaNguard::checkOneParticleGradients({0.0,0.0,0.0}, "unit test"));
}
}
