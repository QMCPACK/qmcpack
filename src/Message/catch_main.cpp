//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "DeviceManager.h"

#ifdef CATCH_MAIN_HAVE_MPI
#include "Message/Communicate.h"
#endif
#include "Host/OutputManager.h"
#include "MemoryUsage.h"
#include "io/hdf/hdf_error_suppression.h"

// Replacement unit test main function to ensure that MPI is finalized once
// (and only once) at the end of the unit test.

// turn on QMCPACK printout
bool turn_on_output = false;

// AFQMC specific unit test arguments.
std::string UTEST_HAMIL, UTEST_WFN;

int main(int argc, char* argv[])
{
  // Suppress HDF5 warning and error messages.
  qmcplusplus::hdf_error_suppression hide_hdf_errors;
  Catch::Session session;
  using namespace Catch::clara;
  // Build command line parser.
  auto cli = session.cli() |
      Opt(UTEST_HAMIL, "UTEST_HAMIL")["--hamil"]("Hamiltonian file to be used by unit test if applicable.") |
      Opt(UTEST_WFN, "UTEST_WFN")["--wfn"]("Wavefunction file to be used by unit test if applicable.") |
      Opt(turn_on_output)["--turn-on-printout"]("Turn on QMCPACK output manager printout");
  session.cli(cli);
  // Parse arguments.
  int parser_err = session.applyCommandLine(argc, argv);
#ifdef CATCH_MAIN_HAVE_MPI
  mpi3::environment env(argc, argv);
  OHMMS::Controller = new Communicate(env.world());
  if (OHMMS::Controller->rank())
    outputManager.shutOff();
  Communicate node_comm{OHMMS::Controller->NodeComm()};
  // assign accelerators within a node
  qmcplusplus::DeviceManager::initializeGlobalDeviceManager(node_comm.rank(), node_comm.size());
#else
  qmcplusplus::DeviceManager::initializeGlobalDeviceManager(0, 1);
#endif
  if (!turn_on_output)
  {
    qmcplusplus::app_log() << "QMCPACK printout is suppressed. Use --turn-on-printout to see all the printout."
                           << std::endl;
    outputManager.shutOff();
  }
  qmcplusplus::print_mem("Before running tests", qmcplusplus::app_log());
  // Run the tests.
  int result = session.run(argc, argv);
#ifdef CATCH_MAIN_HAVE_MPI
  OHMMS::Controller->finalize();
#endif
  if (parser_err != 0)
  {
    return parser_err;
  }
  else
  {
    return result;
  }
}
