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
#include <catch2/matchers/catch_matchers_string.hpp>
#include <vector>
#include <string>
#include "QMCApp/CommandLineOptions.h"
#include "QMCApp/setupQMCState.h"
#include "Utilities/qmc_common.h"
#include "Message/Communicate.h"
#include "Message/UniformCommunicateError.h"

namespace qmcplusplus
{

TEST_CASE("setupQMCState", "[qmcapp]")
{
  CommandLineOptions options;
  QMCState qmc_state;
  Communicate& comm = *OHMMS::Controller;

  SECTION("default options")
  {
    const char* argv[] = {"qmcpack", "test.xml"};
    options.parse(2, const_cast<char**>(argv));
    setupQMCState(options, comm, qmc_state);
    CHECK_FALSE(qmc_state.dryrun);
    CHECK(qmc_state.io_node == (comm.rank() == 0));
  }

  SECTION("--dryrun and --noprint")
  {
    const char* argv[] = {"qmcpack", "--dryrun", "--noprint", "test.xml"};
    options.parse(4, const_cast<char**>(argv));
    setupQMCState(options, comm, qmc_state);
    CHECK(qmc_state.dryrun);
    CHECK_FALSE(qmc_state.io_node);
  }

  SECTION("--vnode-ranks valid")
  {
    const char* argv[] = {"qmcpack", "--vnode-ranks=1", "test.xml"};
    options.parse(3, const_cast<char**>(argv));
    CHECK_NOTHROW(setupQMCState(options, comm, qmc_state));
  }

  SECTION("--vnode-ranks invalid factor")
  {
    const char* argv[] = {"qmcpack", "--vnode-ranks=99999", "test.xml"};
    options.parse(3, const_cast<char**>(argv));
    CHECK_THROWS_WITH(setupQMCState(options, comm, qmc_state),
                      Catch::Matchers::ContainsSubstring("can only be a factor of the uniform node rank count"));
  }

  SECTION("ensemble cannot accommodate all calculations")
  {
    std::vector<std::string> arg_strings = {"qmcpack"};
    for (int k = 0; k < comm.size() + 2; ++k)
      arg_strings.push_back("test" + std::to_string(k) + ".xml");
    std::vector<char*> argv;
    for (auto& s : arg_strings)
      argv.push_back(s.data());
    options.parse(argv.size(), argv.data());
    CHECK_THROWS_WITH(setupQMCState(options, comm, qmc_state),
                      Catch::Matchers::ContainsSubstring("cannot accommodate all the"));
  }
}

} // namespace qmcplusplus
