#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Utilities/SimpleParser.h"
#include <stdio.h>
#include <string>
#include <sstream>
#include <vector>

using std::string;
using std::vector;
namespace qmcplusplus
{
TEST_CASE("parsewords_empty", "[utilities]")
{
  string input = "";
  vector<string> outlist;
  unsigned int num = parsewords(input.c_str(), outlist);

  REQUIRE(num == 0);
  REQUIRE(outlist.size() == 0);
}

struct ParseCase
{
  string input;          // input string
  vector<string> output; // expected tokens
  string extra_split;    // extra characters to split on

  ParseCase(const string& in) : input(in) {}
  ParseCase(const string& in, const vector<string>& out) : input(in), output(out) {}
  ParseCase(const string& in, const vector<string>& out, const string& extra)
      : input(in), output(out), extra_split(extra)
  {}
};
using ParseCaseVector_t = vector<ParseCase>;


TEST_CASE("parsewords", "[utilities]")
{
  ParseCaseVector_t tlist = {// Input string, list of expected tokens, extra  split characters
                             {
                                 "a",
                                 {"a"},
                             },
                             {
                                 "b=",
                                 {"b"},
                             },
                             {
                                 "=c",
                                 {"c"},
                             },
                             {
                                 "d=e",
                                 {"d", "e"},
                             },
                             {
                                 "f,g",
                                 {"f", "g"},
                             },
                             {
                                 "h\ti,j",
                                 {"h", "i", "j"},
                             },
                             {
                                 "k|m",
                                 {"k|m"},
                             },
                             {"n|o", {"n", "o"}, "|"}};
  for (auto& tc : tlist)
  {
    SECTION(string("Parsing string: ") + tc.input)
    {
      vector<string> outlist;
      unsigned int num = parsewords(tc.input.c_str(), outlist, tc.extra_split);
      REQUIRE(num == tc.output.size());
      REQUIRE(outlist.size() == tc.output.size());
      for (int i = 0; i < tc.output.size(); i++)
      {
        REQUIRE(outlist[i] == tc.output[i]);
      }
    }
  }
}


TEST_CASE("readLine", "[utilities]")
{
  ParseCaseVector_t tlist = {
      // Input string, list of expected tokens, extra  split characters
      {"", {""}},
      {"one", {"one"}},
      {"one\ntwo", {"one", "two"}},
      {"one;two", {"one", "two"}},
      {"one\\\ntwo", {"onetwo"}},
      {"one\\ two", {"one\\ two"}},
      {"12345678901234567890extra", {"1234567890123456789"}}, // assuming bufLen=20 below
  };


  for (auto& tc : tlist)
  {
    SECTION(string("Parsing string: ") + tc.input)
    {
      const int bufLen = 20;
      char buf[bufLen];
      std::istringstream input(tc.input);
      for (int i = 0; i < tc.output.size(); i++)
      {
        char* out = readLine(buf, bufLen, input);
        REQUIRE(buf == tc.output[i]);
        if (i == tc.output.size() - 1)
        {
          REQUIRE(out == NULL);
        }
        else
        {
          REQUIRE(out != NULL);
        }
      }
    }
  }
}


TEST_CASE("getwords", "[utilities]")
{
  ParseCaseVector_t tlist = {
      // Input string, list of expected tokens, extra  split characters
      {"one\n", {"one"}},
      {"one,two\n", {"one", "two"}},
      {"one|two\n", {"one|two"}},
      {"a|b\n", {"a", "b"}, "|"},
  };

  for (auto& tc : tlist)
  {
    SECTION(string("Parsing input: ") + tc.input)
    {
      vector<string> outlist;
      std::istringstream input(tc.input);
      int num = getwords(outlist, input, 0, tc.extra_split);
      REQUIRE(num == tc.output.size());
      REQUIRE(outlist.size() == tc.output.size());
      for (int i = 0; i < tc.output.size(); i++)
      {
        REQUIRE(outlist[i] == tc.output[i]);
      }
    }
  }
}

TEST_CASE("getwordsWithMergedNumbers", "[utilities]")
{
  // 1.0-2.0 -> 1.0 -2.0
  // 105 C 5 Z 0.000000 0.000000 -0.567800-103.884284 -2.142253

  ParseCaseVector_t tlist = {// Input string, list of expected tokens, extra  split characters
                             {"1\n", {"1"}},
                             {"1 2\n", {"1", "2"}},
                             {"1.0 -2.0\n", {"1.0", "-2.0"}},
                             {"1.0-2.0\n", {"1.0", "-2.0"}},
                             {"-1.0-2.0-3.0\n", {"-1.0", "-2.0", "-3.0"}},
                             {"105 C 5 Z 0.000000 0.000000 -0.567800-103.884284 -2.142253\n",
                              {"105", "C", "5", "Z", "0.000000", "0.000000", "-0.567800", "-103.884284", "-2.142253"}}};

  for (auto& tc : tlist)
  {
    SECTION(string("Parsing input: ") + tc.input)
    {
      vector<string> outlist;
      std::istringstream input(tc.input);
      int num = getwordsWithMergedNumbers(outlist, input);
      REQUIRE(num == tc.output.size());
      REQUIRE(outlist.size() == tc.output.size());
      for (int i = 0; i < tc.output.size(); i++)
      {
        REQUIRE(outlist[i] == tc.output[i]);
      }
    }
  }
}

} // namespace qmcplusplus
