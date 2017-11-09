//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, mdewing@anl.gov Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include <iostream>
#include "Utilities/UtilityFunctions.h"
#include <stdio.h>
#include <string>
#include <vector>

using std::string;
using std::vector;
namespace qmcplusplus {

void print_vector(vector<int> &out)
{
  for (int i = 0; i < out.size(); i++) {
    std::cout << out[i] << " ";
  }
  std::cout << std::endl;
}

TEST_CASE("FairDivideLow_one", "[utilities]")
{
  std::vector<int> out;
  FairDivideLow(1,1,out);
  REQUIRE(out.size() == 2);
  REQUIRE(out[0] == 0);
  REQUIRE(out[1] == 1);
}

TEST_CASE("FairDivideLow_two", "[utilities]")
{
  std::vector<int> out;
  FairDivideLow(2,1,out);
  REQUIRE(out.size() == 2);
  REQUIRE(out[0] == 0);
  REQUIRE(out[1] == 2);

  FairDivideLow(2,2,out);
  REQUIRE(out.size() == 3);
  REQUIRE(out[0] == 0);
  REQUIRE(out[1] == 1);
  REQUIRE(out[2] == 2);
}

TEST_CASE("FairDivideLow_three", "[utilities]")
{
  std::vector<int> out;
  FairDivideLow(3,1,out);
  REQUIRE(out.size() == 2);
  REQUIRE(out[0] == 0);
  REQUIRE(out[1] == 3);

  FairDivideLow(3,2,out);
  REQUIRE(out.size() == 3);
  REQUIRE(out[0] == 0);
  REQUIRE(out[1] == 1);
  REQUIRE(out[2] == 3);

  FairDivideLow(3,3,out);
  REQUIRE(out.size() == 4);
  REQUIRE(out[0] == 0);
  REQUIRE(out[1] == 1);
  REQUIRE(out[2] == 2);
  REQUIRE(out[3] == 3);
}

TEST_CASE("FairDivideLow_four", "[utilities]")
{
  std::vector<int> out;
  FairDivideLow(4,2,out);
  REQUIRE(out.size() == 3);
  REQUIRE(out[0] == 0);
  REQUIRE(out[1] == 2);
  REQUIRE(out[2] == 4);
}

}
