//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsPETE/TinyVector.h"

#include <stdio.h>
#include <string>

using std::string;


bool pass = true;
namespace qmcplusplus
{

template<unsigned int D>
void test_tiny_vector()
{
  typedef TinyVector<double, D> vec_t;

  vec_t v1;
  // default constructor sets elements to zero
  for (int i = 0; i < D; i++) {
    REQUIRE(v1[i] == Approx(0.0));
  }

  vec_t v2(1.0);
  // single constructor sets all the elements to that value
  for (int i = 0; i < D; i++) {
    REQUIRE(v2[i] == Approx(1.0));
  }

  // TODO: add optional bounds checks to element access methods
  vec_t v4;
  double sum = 0;
  for (int i = 0; i < D; i++) {
    sum += 1.0*i;
    v4[i] = 1.0*i;
  }

  // Dot product
  double dotp = dot(v2, v4);
  REQUIRE(sum == Approx(dotp));

  // Multiply add
  v1 += 2.0*v4;
  REQUIRE(2.0*sum == dot(v2, v1));
}

template <unsigned int D>
void test_tiny_vector_size_two()
{
  typedef TinyVector<double, D> vec_t;
  vec_t v3(1.0, 2.0);
  REQUIRE(v3[0] == Approx(1.0));
  REQUIRE(v3[1] == Approx(2.0));
  // problem: elements past those explicitly set are undefined
  // in this case, vectors with D > 2 will have undefined elements.
}

TEST_CASE("tiny vector", "[OhmmsPETE]")
{
  test_tiny_vector<1>();
  test_tiny_vector<2>();
  test_tiny_vector<3>();
  test_tiny_vector<4>();
  test_tiny_vector_size_two<2>();
  test_tiny_vector_size_two<3>();
  test_tiny_vector_size_two<4>();
}

}
