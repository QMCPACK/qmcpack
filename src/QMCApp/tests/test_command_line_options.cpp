//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>
#include "QMCApp/CommandLineOptions.h"
#include "Message/Communicate.h"
#include "Message/UniformCommunicateError.h"

namespace qmcplusplus
{

TEST_CASE("CommandLineOptions default state", "[qmcapp]")
{
  CommandLineOptions options;
  CHECK_FALSE(options.get_debug());
  CHECK_FALSE(options.get_dryrun());
  CHECK_FALSE(options.get_noprint());
  CHECK_FALSE(options.get_help());
  CHECK_FALSE(options.get_version());
  CHECK_FALSE(options.should_normal_stop());
  CHECK_FALSE(options.get_enable_timers());
  CHECK(options.get_timer_level().empty());
  CHECK_FALSE(options.get_verbosity().has_value());
  CHECK_FALSE(options.get_vnode_ranks().has_value());
  CHECK(options.get_inputs().empty());
}

TEST_CASE("CommandLineOptions input files parsing", "[qmcapp]")
{
  CommandLineOptions options;

  SECTION("single XML input")
  {
    const char* argv[] = {"qmcpack", "test.xml"};
    options.parse(2, const_cast<char**>(argv));
    CHECK_FALSE(options.should_normal_stop());
    CHECK(options.get_inputs() == std::vector<std::string>{"test.xml"});
  }

  SECTION("multiple XML inputs")
  {
    const char* argv[] = {"qmcpack", "job1.xml", "job2.xml"};
    options.parse(3, const_cast<char**>(argv));
    CHECK_FALSE(options.should_normal_stop());
    CHECK(options.get_inputs() == std::vector<std::string>{"job1.xml", "job2.xml"});
  }

  SECTION("XML inputs from text file")
  {
    const std::string list_filename = "test_input_list_" + std::to_string(OHMMS::Controller->rank()) + ".txt";
    {
      std::ofstream fout(list_filename);
      fout << "from_txt_1.xml\n";
      fout << "from_txt_2.xml\n";
    }

    const char* argv[] = {"qmcpack", list_filename.c_str(), "job_cmd.xml"};
    options.parse(3, const_cast<char**>(argv));
    std::filesystem::remove(list_filename);

    const std::vector<std::string> expected{"from_txt_1.xml", "from_txt_2.xml", "job_cmd.xml"};
    CHECK(options.get_inputs() == expected);
  }

  SECTION("no valid input files throws")
  {
    const char* argv1[] = {"qmcpack"};
    CHECK_THROWS_AS(options.parse(1, const_cast<char**>(argv1)), UniformCommunicateError);

    const char* argv2[] = {"qmcpack", "-debug"};
    CHECK_THROWS_AS(options.parse(2, const_cast<char**>(argv2)), UniformCommunicateError);
  }
}

TEST_CASE("CommandLineOptions boolean options", "[qmcapp]")
{
  CommandLineOptions options;

  SECTION("-debug and --debug")
  {
    const char* argv1[] = {"qmcpack", "-debug", "test.xml"};
    options.parse(3, const_cast<char**>(argv1));
    CHECK(options.get_debug());

    const char* argv2[] = {"qmcpack", "--debug", "test.xml"};
    options.parse(3, const_cast<char**>(argv2));
    CHECK(options.get_debug());
  }

  SECTION("--dryrun")
  {
    const char* argv[] = {"qmcpack", "--dryrun", "test.xml"};
    options.parse(3, const_cast<char**>(argv));
    CHECK(options.get_dryrun());
  }

  SECTION("--noprint")
  {
    const char* argv[] = {"qmcpack", "--noprint", "test.xml"};
    options.parse(3, const_cast<char**>(argv));
    CHECK(options.get_noprint());
  }
}

TEST_CASE("CommandLineOptions timer options", "[qmcapp]")
{
  CommandLineOptions options;

  SECTION("-enable-timers without level")
  {
    const char* argv[] = {"qmcpack", "-enable-timers", "test.xml"};
    options.parse(3, const_cast<char**>(argv));
    CHECK(options.get_enable_timers());
    CHECK(options.get_timer_level().empty());
  }

  SECTION("-enable-timers with level")
  {
    const char* argv[] = {"qmcpack", "-enable-timers=timer_level_fine", "test.xml"};
    options.parse(3, const_cast<char**>(argv));
    CHECK(options.get_enable_timers());
    CHECK(options.get_timer_level() == "timer_level_fine");
  }
}

TEST_CASE("CommandLineOptions verbosity option", "[qmcapp]")
{
  CommandLineOptions options;

  SECTION("-verbosity=low")
  {
    const char* argv[] = {"qmcpack", "-verbosity=low", "test.xml"};
    options.parse(3, const_cast<char**>(argv));
    REQUIRE(options.get_verbosity().has_value());
    CHECK(*options.get_verbosity() == Verbosity::LOW);
  }

  SECTION("-verbosity=high")
  {
    const char* argv[] = {"qmcpack", "-verbosity=high", "test.xml"};
    options.parse(3, const_cast<char**>(argv));
    REQUIRE(options.get_verbosity().has_value());
    CHECK(*options.get_verbosity() == Verbosity::HIGH);
  }

  SECTION("-verbosity=debug")
  {
    const char* argv[] = {"qmcpack", "-verbosity=debug", "test.xml"};
    options.parse(3, const_cast<char**>(argv));
    REQUIRE(options.get_verbosity().has_value());
    CHECK(*options.get_verbosity() == Verbosity::DEBUG);
  }

  SECTION("-verbosity invalid")
  {
    const char* argv[] = {"qmcpack", "-verbosity=super_verbose", "test.xml"};
    CHECK_THROWS_MATCHES(options.parse(3, const_cast<char**>(argv)), UniformCommunicateError,
                         Catch::Matchers::Message("Unknown verbosity level: super_verbose"));
  }

  SECTION("--verbosity high (space-separated)")
  {
    const char* argv[] = {"qmcpack", "--verbosity", "high", "test.xml"};
    options.parse(4, const_cast<char**>(argv));
    REQUIRE(options.get_verbosity().has_value());
    CHECK(*options.get_verbosity() == Verbosity::HIGH);
  }

  SECTION("--verbosity missing argument")
  {
    const char* argv[] = {"qmcpack", "--verbosity"};
    CHECK_THROWS_MATCHES(options.parse(2, const_cast<char**>(argv)), UniformCommunicateError,
                         Catch::Matchers::Message("Missing argument for --verbosity"));
  }
}

TEST_CASE("CommandLineOptions vnode-ranks option", "[qmcapp]")
{
  CommandLineOptions options;

  SECTION("--vnode-ranks=<value>")
  {
    const char* argv[] = {"qmcpack", "--vnode-ranks=4", "test.xml"};
    options.parse(3, const_cast<char**>(argv));
    REQUIRE(options.get_vnode_ranks().has_value());
    CHECK(*options.get_vnode_ranks() == 4);
  }

  SECTION("--vnode-ranks <value>")
  {
    const char* argv[] = {"qmcpack", "--vnode-ranks", "8", "test.xml"};
    options.parse(4, const_cast<char**>(argv));
    REQUIRE(options.get_vnode_ranks().has_value());
    CHECK(*options.get_vnode_ranks() == 8);
  }

  SECTION("--vnode-ranks non-positive")
  {
    const char* argv1[] = {"qmcpack", "--vnode-ranks=0", "test.xml"};
    CHECK_THROWS_MATCHES(options.parse(3, const_cast<char**>(argv1)), UniformCommunicateError,
                         Catch::Matchers::Message("--vnode-ranks must be a positive integer."));

    const char* argv2[] = {"qmcpack", "--vnode-ranks=-2", "test.xml"};
    CHECK_THROWS_MATCHES(options.parse(3, const_cast<char**>(argv2)), UniformCommunicateError,
                         Catch::Matchers::Message("--vnode-ranks must be a positive integer."));
  }

  SECTION("--vnode-ranks missing argument")
  {
    const char* argv[] = {"qmcpack", "--vnode-ranks"};
    CHECK_THROWS_MATCHES(options.parse(2, const_cast<char**>(argv)), UniformCommunicateError,
                         Catch::Matchers::Message("Missing argument for --vnode-ranks"));
  }

  SECTION("--vnode-ranks invalid argument")
  {
    const char* argv[] = {"qmcpack", "--vnode-ranks=notanumber", "test.xml"};
    CHECK_THROWS_MATCHES(options.parse(3, const_cast<char**>(argv)), UniformCommunicateError,
                         Catch::Matchers::Message("Invalid argument for --vnode-ranks"));
  }
}

TEST_CASE("CommandLineOptions help, version, and removed options", "[qmcapp]")
{
  CommandLineOptions options;

  SECTION("--help")
  {
    const char* argv[] = {"qmcpack", "--help"};
    options.parse(2, const_cast<char**>(argv));
    CHECK(options.get_help());
    CHECK(options.should_normal_stop());
    CHECK(options.get_inputs().empty());
  }

  SECTION("-h")
  {
    const char* argv[] = {"qmcpack", "-h"};
    options.parse(2, const_cast<char**>(argv));
    CHECK(options.get_help());
    CHECK(options.should_normal_stop());
  }

  SECTION("--version")
  {
    const char* argv[] = {"qmcpack", "--version"};
    options.parse(2, const_cast<char**>(argv));
    CHECK(options.get_version());
    CHECK(options.should_normal_stop());
    CHECK(options.get_inputs().empty());
  }

  SECTION("--save_wfs (removed option throws)")
  {
    const char* argv[] = {"qmcpack", "--save_wfs", "test.xml"};
    CHECK_THROWS_WITH(
        options.parse(3, const_cast<char**>(argv)),
        Catch::Matchers::ContainsSubstring(
            "ERROR: command line option --save_wfs has been removed. Use save_coefs input tag as described in the "
            "manual."));
  }

  SECTION("--vacuum (removed option throws)")
  {
    const char* argv[] = {"qmcpack", "--vacuum", "test.xml"};
    CHECK_THROWS_WITH(
        options.parse(3, const_cast<char**>(argv)),
        Catch::Matchers::ContainsSubstring(
            "ERROR: command line option --vacuum has been removed. Use vacuum input tag as described in the manual."));
  }

  SECTION("multiple error conditions grouped together")
  {
    const char* argv[] = {"qmcpack", "--save_wfs", "--vacuum"};
    CHECK_THROWS_WITH(
        options.parse(3, const_cast<char**>(argv)),
        Catch::Matchers::ContainsSubstring("ERROR: command line option --save_wfs has been removed.") &&
            Catch::Matchers::ContainsSubstring("ERROR: command line option --vacuum has been removed.") &&
            Catch::Matchers::ContainsSubstring("No valid input file is given."));
  }

  SECTION("unrecognized option")
  {
    const char* argv[] = {"qmcpack", "--unknown-opt", "test.xml"};
    CHECK_THROWS_WITH(
        options.parse(3, const_cast<char**>(argv)),
        Catch::Matchers::ContainsSubstring("ERROR: command line option --unknown-opt is unrecognized."));
  }

  SECTION("multiple unrecognized options and missing inputs")
  {
    const char* argv[] = {"qmcpack", "--verbose", "high", "--vnode-rank", "4"};
    CHECK_THROWS_WITH(
        options.parse(5, const_cast<char**>(argv)),
        Catch::Matchers::ContainsSubstring("ERROR: command line option --verbose is unrecognized.") &&
            Catch::Matchers::ContainsSubstring("ERROR: command line option --vnode-rank is unrecognized.") &&
            Catch::Matchers::ContainsSubstring("No valid input file is given."));
  }

  SECTION("unrecognized option with --help does not output error")
  {
    const char* argv[] = {"qmcpack", "--help", "--bogus"};
    options.parse(3, const_cast<char**>(argv));
    CHECK(options.get_help());
    CHECK(options.should_normal_stop());
  }

  SECTION("print_version_and_usage output")
  {
    const char* argv[] = {"qmcpack", "--help"};
    options.parse(2, const_cast<char**>(argv));
    std::ostringstream oss;
    options.print_version_and_usage(oss);
    CHECK_THAT(oss.str(), Catch::Matchers::ContainsSubstring("QMCPACK version"));
    CHECK_THAT(oss.str(), Catch::Matchers::ContainsSubstring("Usage: qmcpack [options] <input-files.xml>"));
    CHECK_THAT(oss.str(), Catch::Matchers::ContainsSubstring("--help, -h"));
    CHECK_THAT(oss.str(), Catch::Matchers::ContainsSubstring("--version"));
    CHECK_THAT(oss.str(), Catch::Matchers::ContainsSubstring("--dryrun"));
    CHECK_THAT(oss.str(), Catch::Matchers::ContainsSubstring("--noprint"));
    CHECK_THAT(oss.str(), Catch::Matchers::ContainsSubstring("-debug, --debug"));
    CHECK_THAT(oss.str(), Catch::Matchers::ContainsSubstring("--enable-timers[=level]"));
    CHECK_THAT(oss.str(), Catch::Matchers::ContainsSubstring("--verbosity=level"));
    CHECK_THAT(oss.str(), Catch::Matchers::ContainsSubstring("--vnode-ranks=N"));
  }
}

} // namespace qmcplusplus
