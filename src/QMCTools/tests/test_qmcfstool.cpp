//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include <iostream>
#include <vector>

#include "Configuration.h"
#include "QMCTools/QMCFiniteSize/SkParserBase.h"
#include "QMCTools/QMCFiniteSize/SkParserASCII.h"
#include <string>

namespace qmcplusplus
{

typedef QMCTraits::RealType RealType;
typedef QMCTraits::PosType PosType;

TEST_CASE("parse Sk file", "[tools]")
{
  SkParserBase* skparser = new SkParserASCII();
  std::string filename = "simple_Sk.dat";
  skparser->parse(filename);
  std::vector<RealType> sk    = skparser->get_sk_raw();
  std::vector<RealType> skerr = skparser->get_skerr_raw();
  std::vector<PosType> grid   = skparser->get_grid_raw();

  REQUIRE(grid[0][0] == -0.395047639918);
  REQUIRE(grid[0][1] == -0.395047639918);
  REQUIRE(grid[0][2] == -0.395047639918);
  REQUIRE(sk[0]      == 0.209093522105);
  REQUIRE(skerr[0]    == 0.00386792459124);

  int last = sk.size() - 1;
  REQUIRE(grid[last][0] == 3.55542875926);
  REQUIRE(grid[last][1] == 1.97523819959);
  REQUIRE(grid[last][2] == 3.55542875926);
  REQUIRE(sk[last]      == 0.999721051092);
  REQUIRE(skerr[last]    == 0.000614284608774);

  delete skparser;
}

} // namespace qmcplusplus
