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
namespace qmcplusplus {


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
  string input; // input string
  string extra_split; // extra characters to split on
  vector<string> output; // expected tokens

  ParseCase(const string &in) : input(in) {}
};
typedef vector<ParseCase> ParseCaseVector_t;


TEST_CASE("parsewords", "[utilities]")
{
#if 0
  // C++ 11 code
  vector<ParseCase> tlist =
  {
    {"a", {"a"}},
    {"b=", {"b"}},
    {"=c", {"c"}}
  };
  for (auto &tc : tlist) {
    vector<string> outlist;
    unsigned int num = parsewords(tc.input.c_str(), outlist);
    REQUIRE(num == tc.output.size());
    REQUIRE(outlist.size() == tc.output.size());
    for (int i = 0; i < tc.output.size(); i++) {
        REQUIRE(outlist[i] == tc.output[i]);
    }
  }
#endif
  ParseCaseVector_t tlist;
  ParseCase p1("a");
  p1.output.push_back("a");
  tlist.push_back(p1);

  ParseCase p2("b=");
  p2.output.push_back("b");
  tlist.push_back(p2);

  ParseCase p3("=c");
  p3.output.push_back("c");
  tlist.push_back(p3);

  ParseCase p4("d=e");
  p4.output.push_back("d");
  p4.output.push_back("e");
  tlist.push_back(p4);
  ParseCase p5("f,g");
  p5.output.push_back("f");
  p5.output.push_back("g");
  tlist.push_back(p5);

  ParseCase p6("h\ti,j");
  p6.output.push_back("h");
  p6.output.push_back("i");
  p6.output.push_back("j");
  tlist.push_back(p6);

  ParseCase p7("k|m");
  p7.output.push_back("k|m");
  tlist.push_back(p7);

  ParseCase p8("n|o");
  p8.extra_split = "|";
  p8.output.push_back("n");
  p8.output.push_back("o");
  tlist.push_back(p8);


  for (ParseCaseVector_t::iterator tc = tlist.begin(); tc != tlist.end(); tc++) {
    SECTION(tc->input) {
      vector<string> outlist;
      unsigned int num = parsewords(tc->input.c_str(), outlist, tc->extra_split);
      REQUIRE(num == tc->output.size());
      REQUIRE(outlist.size() == tc->output.size());
      for (int i = 0; i < tc->output.size(); i++) {
          REQUIRE(outlist[i] == tc->output[i]);
      }
    }
  }
}


TEST_CASE("readLine", "[utilities]")
{
  ParseCaseVector_t tlist;
  ParseCase p1("one");
  p1.output.push_back("one");
  tlist.push_back(p1);

  ParseCase p2("one\ntwo");
  p2.output.push_back("one");
  p2.output.push_back("two");
  tlist.push_back(p2);

  ParseCase p3("one;two");
  p3.output.push_back("one");
  p3.output.push_back("two");
  tlist.push_back(p3);

  ParseCase p4("one\\\ntwo");
  p4.output.push_back("onetwo");
  tlist.push_back(p4);

  ParseCase p5("");
  p5.output.push_back("");
  tlist.push_back(p5);

  ParseCase p6("one\\ two");
  p6.output.push_back("one\\ two");
  tlist.push_back(p6);

  ParseCase p7("12345678901234567890extra");
  p7.output.push_back("1234567890123456789"); // assuming bufLen = 20 below
  tlist.push_back(p7);


  for (ParseCaseVector_t::iterator tc = tlist.begin(); tc != tlist.end(); tc++) {
    SECTION(tc->input) {
      const int bufLen = 20;
      char buf[bufLen];
      std::istringstream input(tc->input);
      for (int i = 0; i < tc->output.size(); i++) {
        char *out = readLine(buf, bufLen, input);
        REQUIRE(buf == tc->output[i]);
        if (i == tc->output.size()-1) {
            REQUIRE(out == NULL);
        } else {
            REQUIRE(out != NULL);
        }
      }
    }
  }
}


TEST_CASE("getwords", "[utilities]")
{
  ParseCaseVector_t tlist;
  ParseCase p1("one\n");
  p1.output.push_back("one");
  tlist.push_back(p1);

  ParseCase p2("one,two\n");
  p2.output.push_back("one");
  p2.output.push_back("two");
  tlist.push_back(p2);

  ParseCase p3("one|two\n");
  p3.output.push_back("one|two");
  tlist.push_back(p3);

  ParseCase p4("a|b\n");
  p4.output.push_back("a");
  p4.output.push_back("b");
  p4.extra_split = "|";
  tlist.push_back(p4);

  for (ParseCaseVector_t::iterator tc = tlist.begin(); tc != tlist.end(); tc++) {
    SECTION(tc->input) {
      vector<string> outlist;
      std::istringstream input(tc->input);
      int num = getwords(outlist, input, 0, tc->extra_split);
      REQUIRE(num == tc->output.size());
      REQUIRE(outlist.size() == tc->output.size());
      for (int i = 0; i < tc->output.size(); i++) {
        REQUIRE(outlist[i] == tc->output[i]);
      }
    }
  }
}

}
