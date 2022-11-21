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
#include "string_utils.h"

/** \file
 */
namespace qmcplusplus
{

TEST_CASE("string_utils_streamToVec", "[utilities]")
{
  { // empty string
    std::string test_string("");
    auto str_vec = convertStrToVec<std::string>(test_string);
    REQUIRE(str_vec.size() == 0);
  }

  { // whitespace string
    std::string test_string("  ");
    auto str_vec = convertStrToVec<std::string>(test_string);
    REQUIRE(str_vec.size() == 0);
  }

  {
    std::string test_string("abc def 123");
    auto str_vec = convertStrToVec<std::string>(test_string);
    REQUIRE(str_vec.size() == 3);
    CHECK(str_vec[0] == "abc");
    CHECK(str_vec[1] == "def");
    CHECK(str_vec[2] == "123");

    // won't work with int
    CHECK_THROWS_WITH(convertStrToVec<int>(test_string),
                      Catch::Matchers::Contains("Error parsing string 'abc def 123' for type (type_info::name) i."));
  }
}
} // namespace qmcplusplus
