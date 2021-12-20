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
#include "ModernStringUtils.hpp"

/** \file
 */
namespace qmcplusplus
{

TEST_CASE("ModernStringUtils_strToLower", "[utilities]")
{
  std::string test_string("IamAMixedCaseTest_String");
  std::string lcase_string = lowerCase(test_string);
  CHECK(lcase_string == "iamamixedcasetest_string");
  std::string i_am_not_just_ascii{"\xab"
                                  "not_Just_ASCII"};
  // very imcomplete check if we are in a c locale that clobbers char beyond ASII
  // \todo If you know C locales well is this ever a worry?
  i_am_not_just_ascii = lowerCase(i_am_not_just_ascii);
  CHECK(i_am_not_just_ascii ==
        "\xab"
        "not_just_ascii");
}
} // namespace qmcplusplus
