//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"
#include <cstdio>
#include <vector>
#include <string>

#include "OhmmsPETE/OhmmsArray.h"

using std::string;

namespace qmcplusplus
{
TEST_CASE("array", "[OhmmsPETE]")
{
  using Array1D = Array<double, 1>;
  Array1D A(3);
  Array1D B(3);

  // iterator
  auto ia = A.begin();
  for (; ia != A.end(); ia++)
  {
    *ia = 1.0;
  }

  // Assignment.   This eventually generates a call to 'evaluate' in OhmmsVector.h
  //  To do: pointer to tutorial on expression template techniques
  B = A;

  for (auto& element : B)
    element *= 3.1;

  CHECK(B(0) == Approx(3.1));
  CHECK(B(1) == Approx(3.1));
  CHECK(B(2) == Approx(3.1));
  REQUIRE(B == B);
  REQUIRE(!(B == A));
  REQUIRE(A != B);
}

TEST_CASE("array NestedContainers", "[OhmmsPETE]")
{
  Array<std::vector<int>, 1> vec_of_vecs({2});
  vec_of_vecs(0).push_back(123);
  vec_of_vecs.resize(5);
  vec_of_vecs(0).clear();
  vec_of_vecs(0).push_back(123);
  vec_of_vecs.resize(0);
  vec_of_vecs.resize(3);
  vec_of_vecs(0).push_back(123);
  CHECK(vec_of_vecs(0).back() == 123);

  Array<std::vector<int>, 1> vec_copy(vec_of_vecs);
  REQUIRE(vec_copy.size() == 3);
  REQUIRE(vec_copy(0).size() == 1);
  CHECK(vec_copy(0).back() == 123);

  Array<std::vector<int>, 1> vec_assign;
  vec_assign = vec_of_vecs;
  REQUIRE(vec_copy.size() == 3);
  REQUIRE(vec_copy(0).size() == 1);
  CHECK(vec_copy(0).back() == 123);
}

TEST_CASE("Array::data", "[OhmmsPETE]")
{
  Array<float, 3> tensor(2, 4, 5);
  REQUIRE(tensor.size() == 40);

  CHECK(tensor.data() + 1 * 4 * 5 + 2 * 5 + 3 == tensor.data_at(1, 2, 3));

  tensor(1, 2, 3) = 0.5f;
  CHECK(*tensor.data_at(1, 2, 3) == 0.5f);
}

TEST_CASE("Array::dimension sizes constructor", "[OhmmsPETE]")
{
  const int dim = 2;
  Array<double, 1> vec(dim);

  Array<double, 3> rank3_tensor(2, 4, 5);
  CHECK(rank3_tensor.shape() == std::array<std::size_t, 3>{2, 4, 5});
  //  rank3_tensor.resize(5,6); this is caught at compile time.
}
} // namespace qmcplusplus
