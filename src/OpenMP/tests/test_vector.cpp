//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <memory>
#include <vector>
#include <iostream>
#include "OpenMP/OMPallocator.hpp"

namespace qmcplusplus
{

TEST_CASE("OMPvector", "[OMP]")
{
  typedef std::vector<double, OMPallocator<double>> vec_t;
  vec_t A(2);

  // iterator
  vec_t::iterator ia = A.begin();
  for (; ia != A.end(); ia++) {
    *ia = 3.1;
  }

  REQUIRE(A[0] == Approx(3.1));
  REQUIRE(A[1] == Approx(3.1));
}

}
