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
#include "Message/CommOperators.h"
#include "Message/Communicate.h"
#include "Utilities/MPIExceptionWrapper.hpp"
#include "Utilities/for_testing/Catch2Approx.h"
#include "VMCBatchedTest.h"

#include <vector>

namespace qmcplusplus
{
TEST_CASE("VMCBatched all-particle moves without drift", "[drivers][mpi]")
{
  MPIExceptionWrapper mpi_exception_wrapper;
  mpi_exception_wrapper([]() {
    Communicate& comm  = *OHMMS::Controller;
    int total_accepted = 0;
    int total_rejected = 0;
    int total_mixed    = 0;
    for (int walkers_per_crowd : {1, 2, 3})
    {
      const auto local = testing::VMCBatchedTest::runAllParticleNoDrift(walkers_per_crowd);
      enum TotalIndex
      {
        WALKERS,
        PROPOSALS,
        RESOLUTIONS,
        ACCEPTED,
        REJECTED,
        MIXED_RESOLUTIONS,
        PARTICLE_MOVES,
        EXPECTED_PARTICLE_MOVES,
        ACCEPTED_PARTICLE_MOVES,
        REJECTED_PARTICLE_MOVES,
        FAILURES,
        TOTAL_COUNT
      };
      std::vector<int> totals(TOTAL_COUNT);
      totals[WALKERS]                 = local.walkers;
      totals[PROPOSALS]               = local.proposals;
      totals[RESOLUTIONS]             = local.resolutions;
      totals[ACCEPTED]                = local.accepted;
      totals[REJECTED]                = local.rejected;
      totals[MIXED_RESOLUTIONS]       = local.mixed_resolutions;
      totals[PARTICLE_MOVES]          = local.particle_moves;
      totals[EXPECTED_PARTICLE_MOVES] = local.expected_particle_moves;
      totals[ACCEPTED_PARTICLE_MOVES] = local.accepted_particle_moves;
      totals[REJECTED_PARTICLE_MOVES] = local.rejected_particle_moves;
      totals[FAILURES]                = static_cast<int>(local.failures.size());
      // Finish collective observations before checking rank-local failures.
      comm.allreduce(totals);
      const int expected_walkers = Concurrency::maxCapacity<>() * walkers_per_crowd * comm.size();
      const int expected_records = testing::VMCBatchedTest::moves_per_walker * expected_walkers;
      CAPTURE(comm.rank(), walkers_per_crowd, local.failures);
      CHECK(local.failures.empty());
      CHECK(local.walkers == Concurrency::maxCapacity<>() * walkers_per_crowd);
      CHECK(local.proposals == testing::VMCBatchedTest::moves_per_walker * local.walkers);
      CHECK(local.resolutions == testing::VMCBatchedTest::moves_per_walker * local.walkers);
      CHECK(local.accepted + local.rejected == local.resolutions);
      CHECK(local.particle_moves == local.expected_particle_moves);
      CHECK(local.accepted_particle_moves + local.rejected_particle_moves == local.expected_particle_moves);
      CHECK(totals[WALKERS] == expected_walkers);
      CHECK(totals[PROPOSALS] == expected_records);
      CHECK(totals[RESOLUTIONS] == expected_records);
      CHECK(totals[ACCEPTED] + totals[REJECTED] == totals[RESOLUTIONS]);
      CHECK(totals[PARTICLE_MOVES] == totals[EXPECTED_PARTICLE_MOVES]);
      CHECK(totals[ACCEPTED_PARTICLE_MOVES] + totals[REJECTED_PARTICLE_MOVES] == totals[EXPECTED_PARTICLE_MOVES]);
      CHECK(totals[FAILURES] == 0);
      // Acceptance is replicated; block weight is root-only, so checking both on root is sufficient.
      if (comm.rank() == 0)
      {
        const double expected_ratio = static_cast<double>(totals[ACCEPTED_PARTICLE_MOVES]) /
            (totals[ACCEPTED_PARTICLE_MOVES] + totals[REJECTED_PARTICLE_MOVES]);
        const double expected_weight = testing::VMCBatchedTest::steps * expected_walkers;
        CHECK(local.driver_acceptance_ratio == Catch::Approx(expected_ratio));
        CHECK(local.driver_block_weight == Catch::Approx(expected_weight));
      }
      total_accepted += totals[ACCEPTED];
      total_rejected += totals[REJECTED];
      total_mixed += totals[MIXED_RESOLUTIONS];
    }
    CHECK(total_accepted > 0);
    CHECK(total_rejected > 0);
    CHECK(total_mixed > 0);
  });
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
