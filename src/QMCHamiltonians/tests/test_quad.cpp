//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Configuration.h"
#include "Numerics/Quadrature.h"

namespace qmcplusplus
{

TEST_CASE("CheckSphericalIntegration", "[hamiltonian]")
{
  // Use the built-in quadrature rule check
  for (int quadrature_index = 1; quadrature_index < 8; quadrature_index++)
  {
    Quadrature3D<QMCTraits::RealType> myRule(quadrature_index,false);
    REQUIRE(myRule.quad_ok);
  }
}

}
