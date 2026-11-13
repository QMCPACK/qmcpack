//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
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
#include "VMCBatchedTest.h"

#include <vector>

namespace qmcplusplus
{
TEST_CASE("VMCBatched all-particle observations across MPI ranks", "[drivers][mpi]")
{
  MPIExceptionWrapper mpi_exception_wrapper;
  mpi_exception_wrapper([]() {
    Communicate& comm = *OHMMS::Controller;
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
      FAILURES,
      TOTAL_COUNT
    };
    int total_accepted = 0;
    int total_rejected = 0;
    int total_mixed    = 0;
    for (int walkers_per_crowd = 1; walkers_per_crowd <= 3; ++walkers_per_crowd)
    {
      const auto local = testing::VMCBatchedTest::runAllParticleNoDrift(walkers_per_crowd);
      std::vector<int> totals(TOTAL_COUNT);
      totals[WALKERS]                 = local.walkers;
      totals[PROPOSALS]               = local.proposals;
      totals[RESOLUTIONS]             = local.resolutions;
      totals[ACCEPTED]                = local.accepted;
      totals[REJECTED]                = local.rejected;
      totals[MIXED_RESOLUTIONS]       = local.mixed_resolutions;
      totals[PARTICLE_MOVES]          = local.particle_moves;
      totals[EXPECTED_PARTICLE_MOVES] = local.expected_particle_moves;
      totals[FAILURES]                = static_cast<int>(local.failures.size());
      // Reduce observations before checking rank-local failures.
      comm.allreduce(totals);
      const int expected_walkers = Concurrency::maxCapacity<>() * walkers_per_crowd * comm.size();
      const int expected_records = testing::VMCBatchedTest::moves_per_walker * expected_walkers;
      CAPTURE(comm.rank(), walkers_per_crowd, local.failures);
      CHECK(local.failures.empty());
      CHECK(totals[WALKERS] == expected_walkers);
      CHECK(totals[PROPOSALS] == expected_records);
      CHECK(totals[RESOLUTIONS] == expected_records);
      CHECK(totals[ACCEPTED] + totals[REJECTED] == totals[RESOLUTIONS]);
      CHECK(totals[PARTICLE_MOVES] == totals[EXPECTED_PARTICLE_MOVES]);
      CHECK(totals[FAILURES] == 0);
      total_accepted += totals[ACCEPTED];
      total_rejected += totals[REJECTED];
      total_mixed += totals[MIXED_RESOLUTIONS];
    }
    CHECK(total_accepted > 0);
    CHECK(total_rejected > 0);
    CHECK(total_mixed > 0);
  });
}
} // namespace qmcplusplus
