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

#include "catch.hpp"

#include "Concurrency/UtilityFunctions.hpp"

namespace qmcplusplus
{
using namespace Concurrency;

TEST_CASE("UtilityFunctions OpenMP", "[concurrency]")
{
  const int old_threads = maxCapacity<Executor::OPENMP>();
  {
    OverrideMaxCapacity<Executor::OPENMP> lock(97);
    CHECK(maxCapacity<Executor::OPENMP>() == 97);
    {
      ThreadCountProtector<Executor::OPENMP> protector;
      omp_set_num_threads(17);
      CHECK(maxCapacity<Executor::OPENMP>() == 17);
    }
    CHECK(maxCapacity<Executor::OPENMP>() == 97);

    {
      OverrideMaxCapacity<Executor::OPENMP> lock(37);
      CHECK(maxCapacity<Executor::OPENMP>() == 37);
    }
    CHECK(maxCapacity<Executor::OPENMP>() == 97);
  }
  CHECK(maxCapacity<Executor::OPENMP>() == old_threads);
}

} // namespace qmcplusplus
