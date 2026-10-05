//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//////////////////////////////////////////////////////////////////////////////////////

#include "ParseCommand.h"
#include "ParserClass.h"

#include <catch2/catch_test_macros.hpp>

#include <cstdio>
#include <fstream>
#include <list>
#include <string>

namespace
{
class TemporaryFile
{
public:
  TemporaryFile(const std::string& name, const std::string& contents) : name_(name)
  {
    std::ofstream output(name_);
    output << contents;
  }

  ~TemporaryFile() { std::remove(name_.c_str()); }

  const std::string& name() const { return name_; }

private:
  std::string name_;
};
} // namespace

// Verify every parser operation reports failure instead of reading beyond an empty buffer.
TEST_CASE("MemParser handles empty input", "[ppconvert][parser]")
{
  TemporaryFile input("ppconvert_empty_input.txt", "");
  MemParserClass parser;
  REQUIRE(parser.OpenFile(input.name()));

  std::string text;
  int integer;
  double real;
  CHECK_FALSE(parser.ReadWord(text));
  CHECK_FALSE(parser.ReadLine(text));
  CHECK_FALSE(parser.NextLine());
  CHECK_FALSE(parser.ReadInt(integer));
  CHECK_FALSE(parser.ReadDouble(real));
}

// Verify tokens without trailing whitespace or a newline retain their final character and stop cleanly at EOF.
TEST_CASE("MemParser reads values ending at EOF", "[ppconvert][parser]")
{
  SECTION("word")
  {
    TemporaryFile input("ppconvert_final_word.txt", "final");
    MemParserClass parser;
    REQUIRE(parser.OpenFile(input.name()));

    std::string word;
    CHECK(parser.ReadWord(word));
    CHECK(word == "final");
    CHECK_FALSE(parser.ReadWord(word));
  }

  SECTION("line")
  {
    TemporaryFile input("ppconvert_final_line.txt", "final line");
    MemParserClass parser;
    REQUIRE(parser.OpenFile(input.name()));

    std::string line;
    CHECK(parser.ReadLine(line));
    CHECK(line == "final line");
    CHECK_FALSE(parser.ReadLine(line));
  }

  SECTION("number")
  {
    TemporaryFile input("ppconvert_final_number.txt", "42");
    MemParserClass parser;
    REQUIRE(parser.OpenFile(input.name()));

    int value;
    CHECK(parser.ReadInt(value));
    CHECK(value == 42);
    CHECK_FALSE(parser.ReadInt(value));
  }
}

// Verify unknown long options are rejected rather than silently added to the accepted argument map.
TEST_CASE("CommandLineParser rejects invalid options", "[ppconvert][command-line]")
{
  std::list<ParamClass> parameters{ParamClass("known", false)};
  CommandLineParserClass parser(parameters);
  char program[]        = "ppconvert";
  char invalid_option[] = "--invalid";
  char* arguments[]     = {program, invalid_option};

  CHECK_FALSE(parser.Parse(2, arguments));
}