//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "OhmmsPETE/TinyVector.h"
#include "checkVector.hpp"
#include <iostream>
#include <string>
#include "OhmmsPETE/OhmmsVector.h"

/** \file
 *  Tests checkVector
 */
namespace qmcplusplus
{

TEST_CASE("checkVector_OhmmsVector_real", "[utilities][for_testing]")
{
  Vector<double> a_vec;
  a_vec.resize(3);
  a_vec[0] = 2.3;
  a_vec[1] = 4.5;
  a_vec[2] = 2.6;

  Vector<double> b_vec;
  b_vec.resize(3);
  b_vec[0] = 2.3;
  b_vec[1] = 4.5;
  b_vec[2] = 2.6;

  auto check_vector_result = checkVector(a_vec, b_vec);
  // This would be how you would fail and print the information about what element failed.
  CHECKED_ELSE(check_vector_result.result) { FAIL(check_vector_result.result_message); }

  //check with epsilon
  b_vec[2] = 2.6005;
  check_vector_result = checkVector(a_vec, b_vec, false, 1e-4);
  REQUIRE(check_vector_result.result == false);
  check_vector_result = checkVector(a_vec, b_vec, false, 1e-3);
  REQUIRE(check_vector_result.result == true);

  b_vec.resize(4);
  b_vec[0] = 2.3;
  b_vec[1] = 4.5;
  b_vec[2] = 2.6;

  check_vector_result = checkVector(a_vec, b_vec);
  CHECKED_ELSE(check_vector_result.result) { FAIL(check_vector_result.result_message); }

  b_vec[0] = 1.0;
  b_vec[1] = 3.6;
  check_vector_result = checkVector(a_vec, b_vec, true);
  std::string::size_type pos = 0;
  int lines                  = 0;
  while (true)
  {
    pos = check_vector_result.result_message.find("\n", pos);
    if (pos == std::string::npos)
      break;
    ++pos;
    ++lines;
  }
  CHECK(lines == 2);
  REQUIRE(check_vector_result.result == false);
}

} // namespace qmcplusplus
