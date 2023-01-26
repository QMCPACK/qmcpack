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


#include "catch.hpp"

#include <memory>
#include <vector>
#include <iostream>
#include "OMPTarget/OMPallocator.hpp"

namespace qmcplusplus
{
TEST_CASE("OMPmath", "[OMP]")
{
  using vec_t = std::vector<double, OMPallocator<double>>;
  vec_t A(3);

  // iterator
  vec_t::iterator ia = A.begin();
  for (; ia != A.end(); ia++)
  {
    *ia = 3.1;
  }

  auto* A_ptr = A.data();
  PRAGMA_OFFLOAD("omp target teams distribute map(always, tofrom:A_ptr[0:2])")
  for (int i = 0; i < 2; i++)
  {
    float s, c, v = 1.2;
    s = std::sin(i * v);
    c = std::cos(i * v);
    //sincos(i*v, &s, &c);
    A_ptr[i] += s + c;
  }

  CHECK(A[0] == Approx(4.1));
  CHECK(A[1] == Approx(4.3943968404));
}

} // namespace qmcplusplus
