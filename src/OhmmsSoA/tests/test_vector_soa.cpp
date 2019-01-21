//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <OhmmsSoA/Container.h>

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{


TEST_CASE("vector", "[OhmmsSoA]")
{
  typedef ParticleAttrib<TinyVector<double,3>> vec_t;
  typedef VectorSoaContainer<double,3> vec_soa_t;

  vec_t R(4);
  vec_soa_t RSoA(4);

  R[0] = TinyVector<double,3>(0.00000000, 0.00000000, 0.00000000);
  R[1] = TinyVector<double,3>(1.68658058, 1.68658058, 1.68658058);
  R[2] = TinyVector<double,3>(3.37316115, 3.37316115, 0.00000000);
  R[3] = TinyVector<double,3>(5.05974172, 5.05974172, 1.68658058);

  RSoA.copyIn(R);

  REQUIRE(RSoA[1][0] == Approx(1.68658058));
  REQUIRE(RSoA[1][1] == Approx(1.68658058));
  REQUIRE(RSoA[1][2] == Approx(1.68658058));

}

}
