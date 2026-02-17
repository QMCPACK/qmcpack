//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include <functional>
#include "catch.hpp"

#include "test_WalkerControl.h"
#include "Message/Communicate.h"
#include "ParticleSetPool.h"
#include "WaveFunctionPool.h"
#include "HamiltonianPool.h"
#include "QMCDrivers/MCPopulation.h"
#include "Utilities/MPIExceptionWrapper.hpp"
#include "Utilities/for_testing/NativeInitializerPrint.hpp"
#include "Platforms/Host/OutputManager.h"


//#include "Concurrency/Info.hpp"
//#include "Concurrency/UtilityFunctions.hpp"

namespace qmcplusplus
{
namespace testing
{
UnifiedDriverWalkerControlMPITest::UnifiedDriverWalkerControlMPITest() : wc_(dpools_.comm, Random)
{
  int num_ranks = dpools_.comm->size();
  pop_ =
      std::make_unique<MCPopulation>(num_ranks, dpools_.comm->rank(), dpools_.particle_pool->getParticleSet("e"),
                                     dpools_.wavefunction_pool->getWaveFunction(), dpools_.hamiltonian_pool->getPrimary());

  pop_->createWalkers(1, walker_confs);
}

void UnifiedDriverWalkerControlMPITest::reportWalkersPerRank(Communicate* c, MCPopulation& pop)
{
#if !defined(NDEBUG)
  std::vector<int> rank_walker_count(c->size(), 0);
  rank_walker_count[c->rank()] = pop.get_num_local_walkers();
  c->allreduce(rank_walker_count);

  const int current_population = std::accumulate(rank_walker_count.begin(), rank_walker_count.end(), 0);

  app_log() << "Walkers Per Rank (Total: " << current_population << ")" << std::endl;
  app_log() << "Rank   Count" << std::endl << "===========" << std::endl;
  for (int i = 0; i < rank_walker_count.size(); ++i)
    app_log() << std::setw(4) << i << "   " << rank_walker_count[i] << std::endl;
#endif
}

void UnifiedDriverWalkerControlMPITest::testNewDistribution(std::vector<int>& initial_num_per_rank,
                                                            std::vector<int>& minus,
                                                            std::vector<int>& plus)
{
  int num_ranks = dpools_.comm->size();
  assert(initial_num_per_rank.size() == num_ranks);
  std::vector<int> num_per_rank = initial_num_per_rank;
  std::vector<int> fair_offset;
  wc_.determineNewWalkerPopulation(num_per_rank, fair_offset, minus, plus);
}

} // namespace testing

TEST_CASE("WalkerControl::determineNewWalkerPopulation", "[drivers][walker_control]")
{
  // These must be empty vectors!
  std::vector<int> minus;
  std::vector<int> plus;
  testing::UnifiedDriverWalkerControlMPITest test;
  int num_ranks = test.getNumRanks();
  std::vector num_per_rank(num_ranks, 1);
  num_per_rank[0] = num_ranks;
  test.testNewDistribution(num_per_rank, minus, plus);
  int rank = test.getRank();
  app_log() << "rank:" << rank << " minus: " << NativePrint(minus) << std::endl;

  std::cout << "rank:" << rank << " plus: " << NativePrint(plus) << std::endl;
  CHECK(minus.size() == num_ranks - 1);
  CHECK(plus.size() == num_ranks - 1);
  app_log() << "rank:" << rank << " plus: " << NativePrint(num_per_rank) << std::endl;
}

void testing::UnifiedDriverWalkerControlMPITest::testPopulationDiff(std::vector<int>& rank_counts_before,
                                                                    std::vector<int>& rank_counts_after)
{
  using MCPWalker = MCPopulation::MCPWalker;

  int rank               = dpools_.comm->rank();
  auto num_local_walkers = pop_->get_num_local_walkers();

  // if rank counts_before is different manipulation multiplicities on walkers to get the rank count we want.
  if (num_local_walkers != rank_counts_before[rank])
  {
    auto walker_diff       = rank_counts_before[rank] - num_local_walkers;
    MCPWalker& last_walker = *(pop_->get_walkers()[num_local_walkers - 1]);
    last_walker.Multiplicity += walker_diff;
  }

  std::vector<int> fair_offset;
  std::vector<int> minus;
  std::vector<int> plus;

  wc_.setNumPerRank(rank_counts_before);
  // note this ordering of calls is from the main code
  wc_.swapWalkersSimple(*pop_);
  wc_.killDeadWalkersOnRank(*pop_);
  pop_->fissionHighMultiplicityWalkers();
  reportWalkersPerRank(dpools_.comm, *pop_);
  CHECK(pop_->get_num_local_walkers() == rank_counts_after[rank]);
}

void testing::UnifiedDriverWalkerControlMPITest::testWalkerIDs(std::vector<std::vector<int>> walker_ids_after,
                                                               std::vector<std::vector<int>> parent_ids_after)
{
  int rank = dpools_.comm->rank();
  std::vector<int> walker_ids;
  std::vector<int> parent_ids;
#ifndef NDEBUG
  for (int iw = 0; iw < pop_->get_walkers().size(); ++iw)
  {
    walker_ids.push_back(pop_->get_walkers()[iw]->getWalkerID());
    parent_ids.push_back(pop_->get_walkers()[iw]->getParentID());
  }
  std::cout << "rank: " << rank << "  walker ids: " << NativePrint(walker_ids)
            << " parent ids: " << NativePrint(parent_ids) << std::endl;
#endif
  for (int iw = 0; iw < walker_ids_after[rank].size(); ++iw)
  {
    CHECK(pop_->get_walkers()[iw]->getWalkerID() == walker_ids_after[rank][iw]);
  }
  for (int iw = 0; iw < parent_ids_after[rank].size(); ++iw)
  {
    CHECK(pop_->get_walkers()[iw]->getParentID() == parent_ids_after[rank][iw]);
  }
}

#define CALC_ID_LAMBDA \
  [num_ranks](int rank, int index_walker_created) { return index_walker_created * num_ranks + rank + 1; }

TEST_CASE("MPI WalkerControl population swap walkers", "[drivers][walker_control]")
{
  auto test_func = []() {
    outputManager.pause();
    testing::UnifiedDriverWalkerControlMPITest test;
    outputManager.resume();

    SECTION("Balanced")
    {
      std::vector<int> count_before(test.getNumRanks(), 1);
      std::vector<int> count_after(test.getNumRanks(), 1);
      // One walker on every node, should be no swapping
      test.testPopulationDiff(count_before, count_after);
    }

    SECTION("LoadBalance")
    {
      auto num_ranks = test.getNumRanks();
      std::vector<int> count_before(num_ranks, 1);
      count_before[0] = num_ranks;
      std::vector<int> count_after(num_ranks, 2);
      count_after[0] = 1;
      app_log() << "count_before: " << NativePrint(count_before) << "  count_after: " << NativePrint(count_after)
                << std::endl;
      test.testPopulationDiff(count_before, count_after);
      std::vector<std::vector<int>> ar_wids;
      std::vector<std::vector<int>> ar_pids;
      auto calcID = CALC_ID_LAMBDA;
      for (int ir = 0; ir < num_ranks; ++ir)
      {
        ar_wids.push_back({ir + 1});
        if (ir > 0)
          ar_wids.back().push_back(calcID(ir, 1));

        ar_pids.push_back({0});
        if (ir > 0)
        {
          ar_pids.back().push_back(1);
        }
      }
      test.testWalkerIDs(ar_wids, ar_pids);
    }

    SECTION("LoadBalance Multiple Copy Optimization")
    {
      int multiple   = 3;
      auto num_ranks = test.getNumRanks();
      std::vector<int> walker_multiplicity_total(num_ranks, 1);
      walker_multiplicity_total[0] = num_ranks * multiple;
      int total_walkers = std::accumulate(walker_multiplicity_total.begin(), walker_multiplicity_total.end(), 0);
      std::vector<int> count_after = fairDivide(total_walkers, num_ranks);
      std::reverse(count_after.begin(), count_after.end());
      app_log() << "walker_multiplicity_before: " << NativePrint(walker_multiplicity_total)
                << "count_after: " << NativePrint(count_after) << std::endl;
      test.testPopulationDiff(walker_multiplicity_total, count_after);
      std::vector<std::vector<int>> ar_wids;
      std::vector<std::vector<int>> ar_pids;
      auto calcID = CALC_ID_LAMBDA;

      for (int ir = 0; ir < num_ranks; ++ir)
      {
        ar_wids.push_back({calcID(ir, 0)});
        for (int iw = 1; iw < count_after[ir]; ++iw)
          ar_wids.back().push_back(calcID(ir, iw));

        ar_pids.push_back({0});
        if (ir == 0)
        {
          for (int iw = 1; iw < count_after[ir]; ++iw)
            ar_pids.back().push_back(calcID(0, 0));
        }
        else
        {
          ar_pids.back().push_back(1);
          for (int iw = 2; iw < count_after[ir]; ++iw)
            ar_pids.back().push_back(calcID(ir, 1));
        }
      }

      app_log() << std::endl;
      test.testWalkerIDs(ar_wids, ar_pids);
    }

    SECTION("MultiStage LoadBalance")
    {
      auto num_ranks = test.getNumRanks();
      std::vector<int> walker_multiplicity_total(num_ranks, 1);
      walker_multiplicity_total[0] = num_ranks;
      std::vector<int> count_after(num_ranks, 2);
      count_after[0] = 1;
      test.testPopulationDiff(walker_multiplicity_total, count_after);
      std::vector<std::vector<int>> ar_wids;
      std::vector<std::vector<int>> ar_pids;
      auto calcID = CALC_ID_LAMBDA;
      for (int ir = 0; ir < num_ranks; ++ir)
      {
        ar_wids.push_back({ir + 1});
        if (ir > 0)
          ar_wids.back().push_back(calcID(ir, 1));

        ar_pids.push_back({0});
        if (ir > 0)
        {
          ar_pids.back().push_back(1);
        }
      }
      test.testWalkerIDs(ar_wids, ar_pids);
      if (num_ranks == 1)
      {
        walker_multiplicity_total = count_after;
      }
      else
      {
        // Next Stage
        walker_multiplicity_total = count_after;
        walker_multiplicity_total[num_ranks - 1] += 2;
        count_after[0] += 1;
        count_after.back() += 1;
        app_log() << "walker_multiplicity_total: " << NativePrint(walker_multiplicity_total)
                  << "  count_after: " << NativePrint(count_after) << std::endl;
        test.testPopulationDiff(walker_multiplicity_total, count_after);
        ar_wids.back().push_back(calcID(num_ranks - 1, 2));
        ar_wids[0].push_back(calcID(0, 1));
        ar_pids.back().push_back(calcID(num_ranks - 1, 1));
        ar_pids[0].push_back(calcID(num_ranks - 1, 1));
        test.testWalkerIDs(ar_wids, ar_pids);

        // Walker dropped
        walker_multiplicity_total = count_after;
        walker_multiplicity_total[num_ranks - 1] -= 1;
        count_after[num_ranks - 1] -= 1;
        test.testPopulationDiff(walker_multiplicity_total, count_after);
        ar_wids.back().pop_back();
        ar_pids.back().pop_back();
        app_log() << "walker_multiplicity_total: " << NativePrint(walker_multiplicity_total)
                  << "  count_after: " << NativePrint(count_after) << std::endl;
        test.testWalkerIDs(ar_wids, ar_pids);

        // Walkers added again
        walker_multiplicity_total = count_after;
        walker_multiplicity_total[num_ranks - 1] += 2;
        count_after[num_ranks - 2] += 1;
        count_after.back() += 1;
        app_log() << "walker_multiplicity_total: " << NativePrint(walker_multiplicity_total)
                  << "  count_after: " << NativePrint(count_after) << std::endl;
        test.testPopulationDiff(walker_multiplicity_total, count_after);
        ar_wids.back().push_back(calcID(num_ranks - 1, 3));
        ar_wids[num_ranks - 2].push_back(calcID(num_ranks - 2, 2));
        ar_pids[num_ranks - 1].push_back(calcID(num_ranks - 1, 1));
        ar_pids[num_ranks - 2].push_back(calcID(num_ranks - 1, 1));
        test.testWalkerIDs(ar_wids, ar_pids);
      }
    }
  };
  MPIExceptionWrapper mew;
  mew(test_func);
}

} // namespace qmcplusplus
