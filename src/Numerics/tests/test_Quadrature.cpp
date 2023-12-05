//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File refactored from QMCHamiltonians/tests/test_ecp.cpp
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"
#include "Configuration.h"
#include "Numerics/Quadrature.h"

namespace qmcplusplus
{
TEST_CASE("CheckSphericalIntegration", "[numerics]")
{
  // Use the built-in quadrature rule check
  for (int quadrature_index = 1; quadrature_index < 9; quadrature_index++)
  {
    Quadrature3D<QMCTraits::RealType> myRule(quadrature_index, false);
    REQUIRE(myRule.quad_ok);
  }
}

} // namespace qmcplusplus
