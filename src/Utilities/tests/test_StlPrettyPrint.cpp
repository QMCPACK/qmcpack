//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include <sstream>
#include "StlPrettyPrint.hpp"

namespace qmcplusplus
{
TEST_CASE("StlPrettyPrint", "[utilities]")
{
  std::ostringstream msg;
  std::vector<int> vec;

  msg << vec;
  CHECK(msg.str() == "[]");

  msg.str(std::string());
  vec = {4, 4, 3, 3};
  msg << vec;
  CHECK(msg.str() == "[4(x2), 3(x2)]");

  msg.str(std::string());
  vec = {4, 4, 3};
  msg << vec;
  CHECK(msg.str() == "[4(x2), 3]");

  msg.str(std::string());
  vec = {4, 3, 3};
  msg << vec;
  CHECK(msg.str() == "[4, 3(x2)]");

  msg.str(std::string());
  vec = {4, 3, 2};
  msg << vec;
  CHECK(msg.str() == "[4, 3, 2]");
}
} // namespace qmcplusplus
