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

#include <catch.hpp>

#include "Message/Communicate.h"
#include "QMCDrivers/DMC/DMCDriverInput.h"
#include "QMCDrivers/DMC/DMCBatched.h"
#include "QMCDrivers/tests/ValidQMCInputSections.h"
#include "QMCDrivers/tests/SetupDMCTest.h"
#include "Concurrency/Info.hpp"
#include "Concurrency/UtilityFunctions.hpp"
#include "Utilities/OutputManager.h"

namespace qmcplusplus
{
namespace testing
{
class DMCBatchedTest
{
public:
  DMCBatchedTest() { up_dtest_ = std::make_unique<SetupDMCTest>(1); }

  void testCalcDefaultLocalWalkers()
  {
    using namespace testing;
    SetupDMCTest& dtest_ = get_dtest();

    auto testWRTWalkersPerRank = [&dtest_](int walkers_per_rank) {
      DMCBatched dmc_batched = dtest_();
      dmc_batched.set_walkers_per_rank(walkers_per_rank, "testing");
      if (dtest_.num_crowds < 8)
        dmc_batched.set_num_crowds(Concurrency::maxThreads(), "Insufficient threads available to match test input");
      QMCDriverNew::AdjustedWalkerCounts awc{0, 0, 0, 0};
      awc.walkers_per_rank = walkers_per_rank;
      awc.num_crowds       = dtest_.num_crowds;
      // what we're testing is whether the awc is transformed to a valid set of walker counts
      awc = dmc_batched.calcDefaultLocalWalkers(awc);

      if (walkers_per_rank < dtest_.num_crowds)
      {
        CHECK(awc.walkers_per_crowd == 1);
        CHECK(awc.num_crowds == awc.walkers_per_rank);
        CHECK(awc.walkers_per_rank == awc.num_crowds);
      }
      else if (walkers_per_rank % dtest_.num_crowds)
      {
        CHECK(awc.walkers_per_crowd == walkers_per_rank / dtest_.num_crowds + 1);
        CHECK(awc.walkers_per_rank == awc.walkers_per_crowd * awc.num_crowds);
      }
      else
      {
        CHECK(awc.walkers_per_rank == awc.walkers_per_crowd * awc.num_crowds);
      }
    };
    testWRTWalkersPerRank(7);
    testWRTWalkersPerRank(31);
    testWRTWalkersPerRank(32);
    testWRTWalkersPerRank(33);
  }


  void testDependentObjectsValidAfterPopulationChange()
  {
    using namespace testing;
    SetupDMCTest& dtest = get_dtest();
  }

  SetupDMCTest& get_dtest() { return *up_dtest_; }

private:
  UPtr<SetupDMCTest> up_dtest_;
};
} // namespace testing

TEST_CASE("DMCBatched::calc_default_local_walkers", "[drivers]")
{
  using namespace testing;

  outputManager.pause();
  DMCBatchedTest dbt;
  outputManager.resume();

  // Its a bit tricky where this needs to be put. Your call needs to be a stack fram below it.
  Concurrency::OverrideMaxThreads<> override(8);
  dbt.testCalcDefaultLocalWalkers();
}

TEST_CASE("DMCBatched change of walkers_per_crowd", "[drivers]")
{
  using namespace testing;

  outputManager.pause();
  DMCBatchedTest dbt;
  outputManager.resume();

  Concurrency::OverrideMaxThreads<> override(8);
}

} // namespace qmcplusplus
