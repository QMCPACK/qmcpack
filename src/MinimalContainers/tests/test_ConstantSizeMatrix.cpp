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

#include "MinimalContainers/ConstantSizeMatrix.hpp"

#include "catch.hpp"

namespace qmcplusplus
{

TEST_CASE("ConstantSizeMatrix stl basics", "[containers]")
{
  ConstantSizeMatrix<double> cmat(1, 9, 1, 32, 0.0);
  CHECK(cmat.size() == 9);
  CHECK(cmat.capacity() == 32);

  CHECK_NOTHROW(cmat.resize(1,16));
  CHECK_THROWS(cmat.resize(1,33));
  CHECK_THROWS(cmat.resize(2,9));

}

}

