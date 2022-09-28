//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"
#include <cstdio>
#include <string>

#include "OhmmsPETE/OhmmsMatrix.h"
#include "config.h"

using std::string;

namespace qmcplusplus
{
TEST_CASE("matrix", "[OhmmsPETE]")
{
  using Mat = Matrix<OHMMS_PRECISION>;
  Mat A(3, 3);
  Mat B(3, 3);
  Mat C(3, 3);

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      B(i, j) = (i + j) * 2.1;
    }
  }

  // Assign to all elements
  A = 1.0;

  // Assignment.   This eventually generates a call to 'evaluate' in OhmmsVector.h
  C = A;

  // *= operator in OhmmMatrixOperators.h
  A *= 3.1;

  // iterator
  Mat::iterator ia = A.begin();
  for (; ia != A.end(); ia++)
  {
    CHECK(*ia == Approx(3.1));
  }

  REQUIRE( A == A);
  REQUIRE( A != C);

  // copy constructor and copy method
  Mat D(C);
  CHECK(D.rows() == 3);
  CHECK(D.cols() == 3);

  // swap_rows
  A(0, 0) = 0.0;
  A(0, 1) = 1.0;
  A(1, 0) = 1.0;
  A(1, 1) = 2.0;
  A.swap_rows(0, 1);
  CHECK(A(0, 0) == 1.0);
  CHECK(A(0, 1) == 2.0);
  CHECK(A(1, 0) == 0.0);
  CHECK(A(1, 1) == 1.0);
  // swap_cols
  A.swap_cols(0, 1);
  CHECK(A(0, 0) == 2.0);
  CHECK(A(0, 1) == 1.0);
  CHECK(A(1, 0) == 1.0);
  CHECK(A(1, 1) == 0.0);
}

TEST_CASE("matrix converting assignment", "[OhmmsPETE]")
{
  Matrix<double> mat_A(3,3);
  Matrix<float> mat_B(3,3);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      mat_A(i, j) = (i + j) * 2.1;

  mat_B = mat_A;

  CHECK(mat_B(0,0) == Approx(0));
  CHECK(mat_B(1,1) == Approx(4.2));

  Matrix<float> mat_C(2,2);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      mat_C(i, j) = (i + j) * 2.2;

  mat_A.assignUpperLeft(mat_C);
  CHECK(mat_A(1,0) == Approx(2.2));
  CHECK(mat_A(1,2) == Approx(6.3));

  Matrix<float> mat_D(4,4);
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      mat_D(i, j) = (i + j) * 2.3;

  mat_A.assignUpperLeft(mat_D);
  CHECK(mat_A(1,0) == Approx(2.3));
  CHECK(mat_A(1,2) == Approx(6.9));
}

} // namespace qmcplusplus
