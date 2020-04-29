//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "Utilities/MPIExceptionWrapper.hpp"
#include "Message/Communicate.h"

namespace qmcplusplus
{
/** Openmp generally works but is not guaranteed with std::atomic
 */
void mpiTestFunctionWrapped(Communicate* comm, std::vector<double> &fake_args)
{
  CHECK(fake_args.size() == 4);
  
  if (comm->size() != 3)
    throw std::runtime_error("Bad Rank Count, test_mpi_exception_wrapper can only be run with 3 MPI ranks.");
}

TEST_CASE("MPIExceptionWrapper function case", "[Utilities]")
{
  Communicate* comm = OHMMS::Controller;

  MPIExceptionWrapper mew;
  std::vector<double> test_vec{1, 2, 3, 4};
  mew(mpiTestFunctionWrapped, comm, test_vec);

  // would really like to check for the MPIAbort but that kills the test program.
}

TEST_CASE("MPIExceptionWrapper lambda case", "[Utilities]")
{
  Communicate* comm = OHMMS::Controller;

  MPIExceptionWrapper mew;
  std::vector<double> test_vec{1, 2, 3, 4};
  auto lambdaTestFunction = [](Communicate* comm, std::vector<double>& fake_args) {
    CHECK(fake_args.size() == 4);
    if (comm->size() != 3)
      throw std::runtime_error("Bad Rank Count, test_mpi_exception_wrapper can only be run with 3 MPI ranks.");
  };
  mew(lambdaTestFunction, comm, test_vec);
}

} // namespace qmcplusplus
