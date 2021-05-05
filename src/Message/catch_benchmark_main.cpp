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

/** \file
 *  \brief enables benchmarking functionality in catch, which has compilation cost.
 *  The only significant difference is the inclusion of catch_benchmark_main.hpp
 */
#define CATCH_CONFIG_RUNNER
#include "catch_benchmark_main.hpp"
#include "Message/Communicate.h"

// Replacement unit test main function to ensure that MPI is finalized once
// (and only once) at the end of the unit test.

// AFQMC specific unit test arguments.
std::string UTEST_HAMIL, UTEST_WFN;

int main(int argc, char* argv[])
{
#ifdef HAVE_MPI
  mpi3::environment env(argc, argv);
  OHMMS::Controller->initialize(env);
#endif
  Catch::Session session;
  using namespace Catch::clara;
  // Build command line parser.
  auto cli = session.cli() |
      Opt(UTEST_HAMIL, "UTEST_HAMIL")["--hamil"]("Hamiltonian file to be used by unit test if applicable.") |
      Opt(UTEST_WFN, "UTEST_WFN")["--wfn"]("Wavefunction file to be used by unit test if applicable.");
  session.cli(cli);
  // Parse arguments.
  int parser_err = session.applyCommandLine(argc, argv);
  // Run the tests.
  int result = session.run(argc, argv);
  OHMMS::Controller->finalize();
  if (parser_err != 0)
  {
    return parser_err;
  }
  else
  {
    return result;
  }
}
