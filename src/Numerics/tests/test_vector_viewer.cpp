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

#include "Numerics/VectorViewer.h"


#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{


TEST_CASE("VectorViewer", "[numerics]")
{
  int a[3];
  a[0] = 2;
  a[1] = 4;
  a[2] = -5;
  VectorViewer<int> view_a(a, 3);

  REQUIRE(view_a.size() == 3);
  
  // operator[]
  REQUIRE(view_a[0] == 2);
  REQUIRE(view_a[1] == 4);
  REQUIRE(view_a[2] == -5);

  // operator[] returning a reference 
  view_a[1] = 42;

  REQUIRE(a[0] == 2);
  REQUIRE(a[1] == 42);
  REQUIRE(a[2] == -5);

 // TODO: add optional bounds checking to accesses via operator[]
}

}
