//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/TinyVector.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
TEST_CASE("vector", "[OhmmsPETE]")
{
  typedef Vector<double> vec_t;
  vec_t A(3);
  vec_t B(3);

  // iterator
  vec_t::iterator ia = A.begin();
  for (; ia != A.end(); ia++)
  {
    *ia = 1.0;
  }

  // Assignment.   This eventually generates a call to 'evaluate' in OhmmsVector.h
  //  To do: pointer to tutorial on expression template techniques
  B = A;

  // *= operator in OhmmVectorOperators.h
  B *= 3.1;

  REQUIRE(B[0] == Approx(3.1));
  REQUIRE(B[1] == Approx(3.1));
  REQUIRE(B[2] == Approx(3.1));
  REQUIRE(B == B);
  REQUIRE(!(B == A));

  vec_t C(2);
  REQUIRE(A != B);
  REQUIRE(A != B);
}

TEST_CASE("Vector simple intializer list", "[OhmmsPETE]")
{
  //empty list should work
  Vector<int> vec_int{};
  Vector<double> vec_double{5.0,4.0,3.0,2.0,1.0};
  CHECK(vec_double[0]==Approx(5.0));
  CHECK(vec_double[4]==Approx(1.0));
}

TEST_CASE("Vector nested intializer list", "[OhmmsPETE]")
{
  Vector<TinyVector<double, 3>> vec_tinyd3({{1,2,3},{4,5,6},{7,8,9}});
  CHECK(vec_tinyd3[1][1] == 5);
  CHECK(vec_tinyd3[2][0] == 7);
  Vector<std::pair<int, int>> vec_pair{{1,2},{3,4}};
  CHECK(vec_pair[0].first == 1);
  CHECK(vec_pair[1].second == 4);
}

TEST_CASE("VectorViewer", "[OhmmsPETE]")
{
  int a[3];
  a[0] = 2;
  a[1] = 4;
  a[2] = -5;
  Vector<int> view_a(a, 3);

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

} // namespace qmcplusplus
