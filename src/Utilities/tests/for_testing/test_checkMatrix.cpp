//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "checkMatrix.hpp"
#include "OhmmsPETE/OhmmsMatrix.h"

namespace qmcplusplus
{
TEST_CASE("checkMatrix_real", "[utilities][for_testing]")
{
  Matrix<double> a_mat;
  a_mat.resize(3, 3);
  a_mat(0, 0) = 2.3;
  a_mat(0, 1) = 4.5;
  a_mat(0, 2) = 2.6;
  a_mat(1, 0) = 0.5;
  a_mat(1, 1) = 8.5;
  a_mat(1, 2) = 3.3;
  a_mat(2, 0) = 1.8;
  a_mat(2, 1) = 4.4;
  a_mat(2, 2) = 4.9;

  Matrix<double> b_mat;
  b_mat.resize(3, 3);
  b_mat(0, 0) = 2.3;
  b_mat(0, 1) = 4.5;
  b_mat(0, 2) = 2.6;
  b_mat(1, 0) = 0.5;
  b_mat(1, 1) = 8.5;
  b_mat(1, 2) = 3.3;
  b_mat(2, 0) = 1.8;
  b_mat(2, 1) = 4.4;
  b_mat(2, 2) = 4.9;

  auto check_matrix_result = checkMatrix(a_mat, b_mat);
  // This would be how you would fail and print the information about what element failed.
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }

  b_mat.resize(4, 4);
  b_mat(0, 0) = 2.3;
  b_mat(0, 1) = 4.5;
  b_mat(0, 2) = 2.6;
  b_mat(1, 0) = 0.5;
  b_mat(1, 1) = 8.5;
  b_mat(1, 2) = 3.3;
  b_mat(2, 0) = 1.8;
  b_mat(2, 1) = 4.4;
  b_mat(2, 2) = 4.9;

  check_matrix_result = checkMatrix(a_mat, b_mat);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }

  b_mat(0, 0) = 1.0;
  b_mat(1,1) = 3.6;
  check_matrix_result = checkMatrix(a_mat, b_mat, true);
  std::cout << check_matrix_result.result_message;
  REQUIRE(check_matrix_result.result == false);
}

} // namespace qmcplusplus
