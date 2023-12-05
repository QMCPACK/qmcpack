//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include <iostream>

#include "Numerics/MinimizeOneDim.h"


namespace qmcplusplus
{
using RealType = double;

class MinTest
{
public:
  MinTest(double value = 0.0) : min_value(value) {}
  RealType min_value;
  RealType one_cycle(RealType x) { return (x - min_value) * (x - min_value); }

  void find_bracket(RealType x0)
  {
    auto bracket = bracket_minimum([this](RealType x) -> RealType { return one_cycle(x); }, x0);
    REQUIRE(bracket.success == true);
    RealType xa = bracket.a;
    RealType xb = bracket.b;
    RealType xc = bracket.c;
    //std::cout << " xa = " << xa;
    //std::cout << " xb = " << xb;
    //std::cout << " xc = " << xc;
    //std::cout << std::endl;

    REQUIRE(xa < xb);
    REQUIRE(xb < xc);

    // For a starting point of 1.3
    //CHECK(xa == Approx(-0.0041));
    //CHECK(xb == Approx( 0.03615));
    //CHECK(xc == Approx(-0.04435));


    RealType fa = one_cycle(xa);
    RealType fb = one_cycle(xb);
    RealType fc = one_cycle(xc);

    REQUIRE(fa > fb);
    REQUIRE(fc > fb);
  }

  // ensure the bracket search will find a minimum at the edge of the bound
  void find_bracket_bound(RealType x0, RealType bound)
  {
    auto bracket = bracket_minimum([this](RealType x) -> RealType { return one_cycle(x); }, x0, bound);
    REQUIRE(bracket.success == false);
  }

  void find_min(RealType x0)
  {
    auto bracket = bracket_minimum([this](RealType x) -> RealType { return one_cycle(x); }, x0);
    auto m       = find_minimum([this](RealType x) -> RealType { return one_cycle(x); }, bracket);

    CHECK(m.first == Approx(min_value));
    CHECK(m.second == Approx(0.0));
  }
};

TEST_CASE("bracket minimum", "[numerics]")
{
  MinTest min_test;
  min_test.find_bracket(1.3);
  min_test.find_bracket(-1.3);
  min_test.find_bracket(10.0);


  MinTest min_test2(1.5);
  min_test2.find_bracket(1.3);
  min_test2.find_bracket(-1.3);
  min_test2.find_bracket(10.0);
  min_test2.find_bracket_bound(1.2, 1.4);

  MinTest min_test3(-0.5);
  min_test3.find_bracket(1.3);
  min_test3.find_bracket(-1.3);
  min_test3.find_bracket(10.0);
  min_test3.find_bracket_bound(1.0, 2.0);
}

TEST_CASE("find minimum", "[numerics]")
{
  MinTest min_test;
  min_test.find_min(1.3);
  min_test.find_min(-1.3);
  min_test.find_min(10.0);

  MinTest min_test2(1.5);
  min_test2.find_min(1.3);
  min_test2.find_min(-1.3);
  min_test2.find_min(10.0);

  MinTest min_test3(-0.5);
  min_test3.find_min(1.3);
  min_test3.find_min(-1.3);
  min_test3.find_min(10.0);
}

} // namespace qmcplusplus
