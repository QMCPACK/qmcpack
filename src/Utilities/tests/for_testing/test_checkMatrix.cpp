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
#include "MatrixAccessor.hpp"
#include <iostream>
#include <string>
#include "OhmmsPETE/OhmmsMatrix.h"

/** \file
 *  Tests checkMatrix and MatrixAccessor.
 *  Separate unit tests are likely a waste of effort  since these are complementary
 *   and in the same module.
 */
namespace qmcplusplus
{
// When linking containers does cause asan failure this should return
TEST_CASE("checkMatrix_OhmmsMatrix_real", "[utilities][for_testing]")
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
  std::string::size_type pos = 0;
  int lines                  = 0;
  while (true)
  {
    pos = check_matrix_result.result_message.find("\n", pos);
    if (pos == std::string::npos)
      break;
    ++pos;
    ++lines;
  }
  CHECK(lines == 2);
  REQUIRE(check_matrix_result.result == false);
}

using testing::MatrixAccessor;

TEST_CASE("checkMatrix_real", "[utilities][for_testing]")
{
  // clang-format off
  std::vector<double> a_vec{2.3, 4.5, 2.6,
                            0.5, 8.5, 3.3,
                            1.8, 4.4, 4.9};
  MatrixAccessor<double> a_mat(a_vec.data(), 3, 3);
  std::vector<double> b_vec{2.3, 4.5, 2.6,
                            0.5, 8.5, 3.3,
                            1.8, 4.4, 4.9};
  MatrixAccessor<double> b_mat(b_vec.data(), 3, 3);

  auto check_matrix_result = checkMatrix(a_mat, b_mat);
  // This would be how you would fail and print the information about what element failed.
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }

  std::vector<double> c_vec{2.3, 4.5, 2.6, 0.0,
                            0.5, 8.5, 3.3, 0.0,
                            1.8, 4.4, 4.9, 0,0,
                            0.0, 0.0, 0.0, 0.0};
  MatrixAccessor<double> c_mat(c_vec.data(), 4, 4);
  // clang-format on

  check_matrix_result = checkMatrix(a_mat, c_mat);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }

  b_vec[0]                   = 1.0;
  b_vec[5]                   = 3.6;
  check_matrix_result        = checkMatrix(a_mat, b_mat, true);
  std::string::size_type pos = 0;
  int lines                  = 0;
  while (true)
  {
    pos = check_matrix_result.result_message.find("\n", pos);
    if (pos == std::string::npos)
      break;
    ++pos;
    ++lines;
  }
  CHECK(lines == 2);
  REQUIRE(check_matrix_result.result == false);
}

} // namespace qmcplusplus
