//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsSoA/VectorSoaContainer.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
TEST_CASE("vector", "[OhmmsSoA]")
{
  using vec_t     = Vector<TinyVector<double, 3>>;
  using vec_soa_t = VectorSoaContainer<double, 3>;

  vec_t R(4);
  vec_soa_t RSoA(4);

  R[0] = TinyVector<double, 3>(0.00000000, 0.00000000, 0.00000000);
  R[1] = TinyVector<double, 3>(1.68658058, 1.68658058, 1.68658058);
  R[2] = TinyVector<double, 3>(3.37316115, 3.37316115, 0.00000000);
  R[3] = TinyVector<double, 3>(5.05974172, 5.05974172, 1.68658058);

  RSoA.copyIn(R);

  //check in value
  CHECK(R[1][0] == Approx(1.68658058));
  CHECK(R[1][1] == Approx(1.68658058));
  CHECK(R[1][2] == Approx(1.68658058));

  //check out value
  CHECK(RSoA[1][0] == Approx(1.68658058));
  CHECK(RSoA[1][1] == Approx(1.68658058));
  CHECK(RSoA[1][2] == Approx(1.68658058));

  //check view value
  vec_soa_t RSoA_view(RSoA.data(), RSoA.size(), RSoA.capacity());
  CHECK(RSoA_view[1][0] == Approx(1.68658058));
  CHECK(RSoA_view[1][1] == Approx(1.68658058));
  CHECK(RSoA_view[1][2] == Approx(1.68658058));
}

TEST_CASE("VectorSoaContainer copy constructor", "[OhmmsSoA]")
{
  Vector<TinyVector<double, 3>> R(4);
  VectorSoaContainer<double, 3> RSoA(4);

  R[0] = TinyVector<double, 3>(0.00000000, 0.00000000, 0.00000000);
  R[1] = TinyVector<double, 3>(1.68658058, 1.68658058, 1.68658058);
  R[2] = TinyVector<double, 3>(3.37316115, 3.37316115, 0.00000000);
  R[3] = TinyVector<double, 3>(5.05974172, 5.05974172, 1.68658058);

  RSoA.copyIn(R);

  VectorSoaContainer<double, 3> rsoa_copy(RSoA);

  // more importantly this test shall not leak memory

  //check out value
  CHECK(rsoa_copy[1][0] == Approx(1.68658058));
  CHECK(rsoa_copy[1][1] == Approx(1.68658058));
  CHECK(rsoa_copy[1][2] == Approx(1.68658058));
}

TEST_CASE("VectorSoaContainer move constructor", "[OhmmsSoA]")
{
  Vector<TinyVector<double, 3>> R(4);
  VectorSoaContainer<double, 3> RSoA(4);

  R[0] = TinyVector<double, 3>(0.00000000, 0.00000000, 0.00000000);
  R[1] = TinyVector<double, 3>(1.68658058, 1.68658058, 1.68658058);
  R[2] = TinyVector<double, 3>(3.37316115, 3.37316115, 0.00000000);
  R[3] = TinyVector<double, 3>(5.05974172, 5.05974172, 1.68658058);

  RSoA.copyIn(R);

  VectorSoaContainer<double, 3> rsoa_move(std::move(RSoA));

  // more importantly this test shall not leak memory

  //check out value
  CHECK(rsoa_move[1][0] == Approx(1.68658058));
  CHECK(rsoa_move[1][1] == Approx(1.68658058));
  CHECK(rsoa_move[1][2] == Approx(1.68658058));
}

TEST_CASE("VectorSoaContainer assignment", "[OhmmsSoA]")
{
  Vector<TinyVector<double, 3>> R(4);
  VectorSoaContainer<double, 3> RSoA(4);

  R[0] = TinyVector<double, 3>(0.00000000, 0.00000000, 0.00000000);
  R[1] = TinyVector<double, 3>(1.68658058, 1.68658058, 1.68658058);
  R[2] = TinyVector<double, 3>(3.37316115, 3.37316115, 0.00000000);
  R[3] = TinyVector<double, 3>(5.05974172, 5.05974172, 1.68658058);
  RSoA.copyIn(R);

  VectorSoaContainer<double, 3> rsoa_assign;
  rsoa_assign = RSoA;
  CHECK(rsoa_assign[3][0] == Approx(5.05974172));
  CHECK(rsoa_assign[3][1] == Approx(5.05974172));
  CHECK(rsoa_assign[3][2] == Approx(1.68658058));

  Vector<TinyVector<double, 3>> r_big(5);
  VectorSoaContainer<double, 3> r_soa_big(5);

  r_big[0] = TinyVector<double, 3>(0.00000000, 0.00000000, 0.00000000);
  r_big[1] = TinyVector<double, 3>(1.68658058, 1.68658058, 1.68658058);
  r_big[2] = TinyVector<double, 3>(3.37316115, 3.37316115, 0.00000000);
  r_big[3] = TinyVector<double, 3>(5.05974172, 5.05974172, 1.68658058);
  r_big[4] = TinyVector<double, 3>(3.37316115, 3.37316115, 0.00000000);

  r_soa_big.copyIn(r_big);

  rsoa_assign = r_soa_big;
  CHECK(rsoa_assign[4][0] == Approx(3.37316115));
  CHECK(rsoa_assign[4][2] == Approx(0.00000000));
}

} // namespace qmcplusplus
