//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include <catch.hpp>

#include "Message/Communicate.h"
#include "QMCDrivers/DMC/DMCDriverInput.h"
#include "QMCDrivers/DMC/DMCBatched.h"
#include "QMCDrivers/tests/ValidQMCInputSections.h"
#include "QMCDrivers/tests/SetupDMCTest.h"
#include "Concurrency/Info.hpp"
#include "Concurrency/UtilityFunctions.hpp"

namespace qmcplusplus
{
TEST_CASE("DMCBatched::calc_default_local_walkers", "[drivers]")
{
  using namespace testing;
  Concurrency::OverrideMaxThreads<> override(8);

  SetupDMCTest dtest(1);


  auto testWRTWalkersPerRank = [&dtest](int walkers_per_rank) {
                                 DMCBatched dmc_batched = dtest();
    dmc_batched.set_walkers_per_rank(walkers_per_rank, "testing");
    if (dtest.num_crowds < 8)
      dmc_batched.set_num_crowds(Concurrency::maxThreads(), "Insufficient threads available to match test input");
    DMCBatched::IndexType local_walkers       = dmc_batched.calc_default_local_walkers(walkers_per_rank);
    QMCDriverNew::IndexType walkers_per_crowd = dmc_batched.get_walkers_per_crowd();

    if (walkers_per_rank < dtest.num_crowds)
    {
      CHECK(walkers_per_crowd == 1);
      CHECK(local_walkers == dtest.num_crowds);
      CHECK(dtest.population.get_num_local_walkers() == dtest.num_crowds);
      CHECK(dtest.population.get_num_global_walkers() == dtest.num_crowds * dtest.num_ranks);
    }
    else if (walkers_per_rank % dtest.num_crowds)
    {
      CHECK(walkers_per_crowd == walkers_per_rank / dtest.num_crowds + 1);
      CHECK(local_walkers == walkers_per_crowd * dtest.num_crowds);
      CHECK(dtest.population.get_num_local_walkers() == walkers_per_crowd * dtest.num_crowds);
      CHECK(dtest.population.get_num_global_walkers() == walkers_per_crowd * dtest.num_crowds * dtest.num_ranks);
    }
    else
    {
      CHECK(local_walkers == walkers_per_rank);
      CHECK(dtest.population.get_num_local_walkers() == walkers_per_rank);
      CHECK(dtest.population.get_num_global_walkers() == walkers_per_rank * dtest.num_ranks);
    }
  };
  testWRTWalkersPerRank(7);
  testWRTWalkersPerRank(31);
  testWRTWalkersPerRank(32);
  testWRTWalkersPerRank(33);
}


} // namespace qmcplusplus
