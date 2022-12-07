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

TEST_CASE("ModernStringUtils_split", "[utilities]")
{
  using modernstrutil::split;

  std::string nothing;
  auto no_tokens = split(nothing, " ");
  CHECK(no_tokens.empty());
  no_tokens = split(nothing, "");
  CHECK(no_tokens.empty());

  std::string white_space{"             "};
  auto tokens = split(white_space, ".");
  CHECK(tokens.size() == 1);
  CHECK(tokens[0] == "             ");

  tokens = split(white_space, " ");
  CHECK(tokens.empty());
  
  std::string test_line{"hi there 101, random line"};
  tokens = split(test_line, " ");
  CHECK(tokens[0].size() == 2);
  CHECK(tokens[4].size() == 4);
  CHECK(tokens[3] == "random");

  tokens = split(test_line, "");
  CHECK(tokens.size() == 1);

  tokens = split(test_line, ";");
  CHECK(tokens.size() == 1);
  
  std::string test_lines{R"(
this is a multi
line
token test
)"};
  auto tokens_lines = split(test_lines, "\n");
  CHECK(tokens_lines[0] == "this is a multi");
  CHECK(tokens_lines[1] == "line");
  CHECK(tokens_lines[2] == "token test");

  std::string test_multidel{"this \t is a multidelimiter  \n   \n token test"};
  auto tokens_multidel = split(test_multidel, "\t \n");
  CHECK(tokens_multidel[0] == "this");
  CHECK(tokens_multidel[1] == "is");
  CHECK(tokens_multidel[4] == "token");
}

TEST_CASE("ModernStringUtils_string2Real", "[utilities]")
{
  std::string_view svalue{"101.52326626"};
  double value = string2Real<double>(svalue);
  CHECK(value == Approx(101.52326626));
} // namespace qmcplusplus

TEST_CASE("ModernStringUtils_strip", "[utilities]")
{
  using modernstrutil::strip;

  std::string nothing;
  auto stripped_nothing = strip(nothing);
  CHECK(stripped_nothing.empty());

  std::string white_space{"           "};
  auto stripped_white_space = strip(white_space);
  CHECK(stripped_white_space.empty());

  std::string test_lines{R"(
    r1 1 0 0
    r2 0 1 0
    r3 0 0 1
)"};

  std::string_view stripped_lines = strip(test_lines);

  std::string_view ref_stripped{"r1 1 0 0\n    r2 0 1 0\n    r3 0 0 1"};
  CHECK(ref_stripped == stripped_lines);

  std::string_view another{"r1 1 0 0\n    r2 0 1 0\n    r3 0 0 1\n  \0"};
  std::string_view another_stripped = strip(another);
  CHECK(another_stripped == "r1 1 0 0\n    r2 0 1 0\n    r3 0 0 1");

  std::string_view unneeded{"this needs no stripping"};
  std::string_view unneeded_stripped = strip(unneeded);
  CHECK(unneeded_stripped == unneeded);
}
  
}
