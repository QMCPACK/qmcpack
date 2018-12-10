//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Numerics/e2iphi.h"


#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{

template <unsigned int N, typename T>
void test_e2iphi()
{
  T phi[N];
  T vcos[N];
  T vsin[N];
  for (int i = 0; i < N; i++) {
    phi[i] = 0.2*i;
  }

  eval_e2iphi(N, phi, vcos, vsin);

  for (int i = 0; i < N; i++) {
    REQUIRE(vcos[i] == Approx(std::cos(phi[i])));
    REQUIRE(vsin[i] == Approx(std::sin(phi[i])));
  }
}

TEST_CASE("e2iphi", "[numerics]")
{
  test_e2iphi<1, double>();
  test_e2iphi<2, double>();
  test_e2iphi<3, double>();
  test_e2iphi<4, double>();

  test_e2iphi<1, float>();
  test_e2iphi<2, float>();
  test_e2iphi<3, float>();
  test_e2iphi<4, float>();
}

}
