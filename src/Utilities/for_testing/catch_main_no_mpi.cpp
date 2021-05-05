//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

int main(int argc, char* argv[])
{
  Catch::Session session;
  using namespace Catch::clara;
  // Build command line parser.
  auto cli = session.cli();
  session.cli(cli);
  // Parse arguments.
  int parser_err = session.applyCommandLine(argc, argv);
  // Run the tests.
  int result = session.run(argc, argv);
  if (parser_err != 0)
  {
    return parser_err;
  }
  else
  {
    return result;
  }
}
