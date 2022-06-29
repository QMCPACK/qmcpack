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

#include "catch.hpp"
#include <iostream>
#include <vector>

#include "Configuration.h"
#include "QMCTools/QMCFiniteSize/SkParserBase.h"
#include "QMCTools/QMCFiniteSize/SkParserASCII.h"
#include "QMCTools/QMCFiniteSize/QMCFiniteSize.h"
#include "Particle/ParticleSet.h"
#include "Particle/ParticleSetPool.h"
#include <string>

namespace qmcplusplus
{
TEST_CASE("FS parse Sk file", "[tools]")
{
  using RealType                         = QMCTraits::RealType;
  using PosType                          = QMCTraits::PosType;
  std::unique_ptr<SkParserBase> skparser = std::make_unique<SkParserASCII>();
  std::string filename                   = "simple_Sk.dat";
  skparser->parse(filename);
  std::vector<RealType> sk    = skparser->get_sk_raw();
  std::vector<RealType> skerr = skparser->get_skerr_raw();
  std::vector<PosType> grid   = skparser->get_grid_raw();

  REQUIRE(grid[0][0] == Approx(0.0));
  REQUIRE(grid[0][1] == Approx(0.0));
  REQUIRE(grid[0][2] == Approx(-6.283185307179586));
  REQUIRE(sk[0] == Approx(0.07225651367144714));
  REQUIRE(skerr[0] == Approx(0.01));

  int last = sk.size() - 1;
  REQUIRE(grid[last][0] == Approx(-18.84955592153876));
  REQUIRE(grid[last][1] == Approx(18.84955592153876));
  REQUIRE(grid[last][2] == Approx(75.39822368615503));
  REQUIRE(sk[last] == Approx(0.9999947116274186));
  REQUIRE(skerr[last] == Approx(0.01));
}

TEST_CASE("FS evaluate", "[tools]")
{
  using RealType = QMCTraits::RealType;
  using PosType  = QMCTraits::PosType;

  std::unique_ptr<SkParserBase> skparser = std::make_unique<SkParserASCII>();
  std::string filename                   = "simple_Sk.dat";
  skparser->parse(filename);

  QMCFiniteSize qfs(skparser.get());
  qfs.parse(std::string("simple_input.xml"));
  qfs.validateXML();
  qfs.initialize();

  /// reference numbers and simple_Sk.dat created by fs_ref.py
  std::vector<RealType> skr = skparser->get_sk_raw();
  RealType vsum             = qfs.calcPotentialDiscrete(skr);
  REQUIRE(vsum == Approx(1.0547517220577185));

  std::vector<RealType> sk, skerr;
  skparser->get_sk(sk, skerr);
  RealType vint = qfs.calcPotentialInt(sk);
  REQUIRE(vint == Approx(1.066688342657357).epsilon(0.001));
}

} // namespace qmcplusplus
