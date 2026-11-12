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
#include <catch2/catch_test_macros.hpp>

#include "Concurrency/Info.hpp"
#include "VMCBatchedTest.h"

namespace qmcplusplus
{
TEST_CASE("VMCBatched all-particle moves without drift", "[drivers]")
{
  int total_accepted = 0;
  int total_rejected = 0;
  int total_mixed    = 0;
  for (int walkers_per_crowd : {1, 2, 3})
  {
    const auto summary         = testing::VMCBatchedTest::runAllParticleNoDrift(walkers_per_crowd);
    const int expected_walkers = Concurrency::maxCapacity<>() * walkers_per_crowd;
    CHECK(summary.walkers == expected_walkers);
    CAPTURE(walkers_per_crowd, summary.failures);
    CHECK(summary.failures.empty());
    CHECK(summary.proposals == testing::VMCBatchedTest::moves_per_walker * expected_walkers);
    CHECK(summary.resolutions == testing::VMCBatchedTest::moves_per_walker * expected_walkers);
    CHECK(summary.accepted + summary.rejected == summary.resolutions);
    CHECK(summary.particle_moves == summary.expected_particle_moves);
    total_accepted += summary.accepted;
    total_rejected += summary.rejected;
    total_mixed += summary.mixed_resolutions;
  }
  CHECK(total_accepted > 0);
  CHECK(total_rejected > 0);
  CHECK(total_mixed > 0);
}

TEST_CASE("VMCBatched all-particle boundary rejection restores the whole walker", "[drivers][vmc]")
{
  for (int invalid_particle : {0, 1})
  {
    CAPTURE(invalid_particle);
    testing::VMCBatchedTest::checkAllParticleBoundaryRejection(invalid_particle);
  }
}

} // namespace qmcplusplus
