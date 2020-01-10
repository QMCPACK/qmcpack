//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "MinimalContainers/ConstantSizeVector.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

namespace qmcplusplus
{

TEST_CASE("ConstantSizeVector stl basics", "[containers]")
{
  ConstantSizeVector<double> cvec(9, 32, 0.0);
  CHECK(cvec.size() == 9);
  CHECK(cvec.capacity() == 32);

  std::vector<double> std_vec(24,1.0);
  cvec = std_vec;
  CHECK(cvec.size() == 24);
  CHECK(cvec.capacity() == 32);
  CHECK(cvec(23)==1.0);

  std_vec.resize(33);
  CHECK_THROWS(cvec = std_vec);
}

}
