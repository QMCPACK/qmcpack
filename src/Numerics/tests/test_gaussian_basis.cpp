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


#include "catch.hpp"
#include "Numerics/GaussianBasisSet.h"

#include <stdio.h>
#include <string>

using std::string;

typedef OHMMS_PRECISION real_type;

namespace qmcplusplus
{

TEST_CASE("Basic Gaussian", "[numerics]") {
  GaussianCombo<real_type>::BasicGaussian g1(3.0, 1.0);

  real_type alpha = 3.0;
  g1.reset(alpha, 1.0);

  real_type r = 1.2;

  real_type f = g1.f(r*r);
  //std::cout << "f = " << f << std::endl << std::endl;
  REQUIRE(f == Approx(0.0132998835424438));

  real_type df = g1.df(r, r*r);
  //std::cout << "df = " << df << std::endl << std::endl;
  REQUIRE(df == Approx(-0.0957591615055951));

  df = 0.0;
  real_type ddf = 0.0;
  f = g1.evaluate(r, r*r, df, ddf);
  //std::cout << "f = " << f << " " << df << " " << ddf << std::endl;
  REQUIRE(f == Approx(0.0132998835424438));
  REQUIRE(df == Approx(-0.0957591615055951));
  REQUIRE(ddf == Approx(0.609666661585622));

  df = 0.0;
  ddf = 0.0;
  real_type d3f = 0.0;
  f = g1.evaluate(r, r*r, df, ddf, d3f);
  //std::cout << "f = " << f << " " << df << " " << ddf  << " " << d3f << std::endl;
  REQUIRE(f == Approx(0.0132998835424438));
  REQUIRE(df == Approx(-0.0957591615055951));
  REQUIRE(ddf == Approx(0.609666661585622));
  REQUIRE(d3f == Approx(-3.24049002534934));
}

TEST_CASE("Gaussian Combo", "[numerics]") {
  GaussianCombo<real_type> gc(0,false);

  // STO-3G for H
  gc.addGaussian(0.15432897, 3.42525091);
  gc.addGaussian(0.53532814, 0.62391373);
  gc.addGaussian(0.44463454, 0.16885540);

  REQUIRE(gc.size() == 3);

  real_type r = 1.3;
  real_type f = gc.f(r);
  //std::cout << "f = " << f << std::endl << std::endl;
  REQUIRE(f == Approx(0.556240444149480));

  f = gc.evaluate(r,1.0/r);
  REQUIRE(f == Approx(0.556240444149480));

  real_type df = gc.df(r);
  REQUIRE(df == Approx(-0.661028435778766));

  gc.evaluateAll(r,1.0/r);
  REQUIRE(gc.Y == Approx(0.556240444149480));
  REQUIRE(gc.dY == Approx(-0.661028435778766));
  REQUIRE(gc.d2Y == Approx(0.643259180749128));

  gc.evaluateWithThirdDeriv(r,1.0/r);
  REQUIRE(gc.Y == Approx(0.556240444149480));
  REQUIRE(gc.dY == Approx(-0.661028435778766));
  REQUIRE(gc.d2Y == Approx(0.643259180749128));
  REQUIRE(gc.d3Y == Approx(-0.896186412781167));
}

TEST_CASE("Gaussian Combo P", "[numerics]") {
  GaussianCombo<real_type> gc(1,false);

  // cc-pVDZ for C, the 2P basis
  gc.addGaussian(0.0381090, 9.4390000);
  gc.addGaussian(0.2094800, 2.0020000);
  gc.addGaussian(0.5085570, 0.5456000);

  REQUIRE(gc.size() == 3);

  real_type r = 1.3;
  real_type f = gc.f(r);

  REQUIRE(f == Approx(0.326057642350121));

  f = gc.evaluate(r,1.0/r);
  REQUIRE(f == Approx(0.326057642350121));

  real_type df = gc.df(r);
  REQUIRE(df == Approx(-0.649531407846947));

  gc.evaluateAll(r,1.0/r);
  REQUIRE(gc.Y == Approx(0.326057642350121));
  REQUIRE(gc.dY == Approx(-0.649531407846947));
  REQUIRE(gc.d2Y == Approx(1.39522444199589));

  gc.evaluateWithThirdDeriv(r,1.0/r);
  REQUIRE(gc.Y == Approx(0.326057642350121));
  REQUIRE(gc.dY == Approx(-0.649531407846947));
  REQUIRE(gc.d2Y == Approx(1.39522444199589));
  REQUIRE(gc.d3Y == Approx(-3.38467690038774));
}

TEST_CASE("Gaussian Combo D", "[numerics]") {
  GaussianCombo<real_type> gc(2,false);

  // cc-pVDZ for C, the 3D basis
  gc.addGaussian(1.0, 0.5500000);

  REQUIRE(gc.size() == 1);

  real_type r = 1.3;
  real_type f = gc.f(r);

  REQUIRE(f == Approx(0.361815669819519));

  f = gc.evaluate(r,1.0/r);
  REQUIRE(f == Approx(0.361815669819519));

  real_type df = gc.df(r);
  REQUIRE(df == Approx(-0.517396407841913));

  gc.evaluateAll(r,1.0/r);
  REQUIRE(gc.Y == Approx(0.361815669819519));
  REQUIRE(gc.dY == Approx(-0.517396407841913));
  REQUIRE(gc.d2Y == Approx(0.341879626412464));

  gc.evaluateWithThirdDeriv(r,1.0/r);
  REQUIRE(gc.Y == Approx(0.361815669819519));
  REQUIRE(gc.dY == Approx(-0.517396407841913));
  REQUIRE(gc.d2Y == Approx(0.341879626412464));
  REQUIRE(gc.d3Y == Approx(0.649384231482385));

}
}
