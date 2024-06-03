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
  // if (num_ranks != 3)
  //   throw std::runtime_error("Bad Rank Count, WalkerControlMPI tests can only be run with 3 MPI ranks.");
  pop_ =
      std::make_unique<MCPopulation>(num_ranks, dpools_.comm->rank(), dpools_.particle_pool->getParticleSet("e"),
                                     dpools_.wavefunction_pool->getPrimary(), dpools_.hamiltonian_pool->getPrimary());

  pop_->createWalkers(1, walker_confs);
}

/** Getting the "fat" walker valid enough to be MPI swapable
 *  
 *  By no means is this "valid" from the perspective of running QMC
 *  See QMCDriverNew::initialLogEvaluation
 *  A fat walker does not seem to be "valid" until all that is done.
 */
void UnifiedDriverWalkerControlMPITest::makeValidWalkers()
{
  auto walker_elements = pop_->get_walker_elements();

  for (auto we : walker_elements)
  {
    we.pset.update();
    if (we.walker.DataSet.size() <= 0)
    {
      we.walker.registerData();
      we.twf.registerData(we.pset, we.walker.DataSet);
      we.walker.DataSet.allocate();
    }
    we.twf.copyFromBuffer(we.pset, we.walker.DataSet);
    we.twf.evaluateLog(we.pset);
    we.twf.updateBuffer(we.pset, we.walker.DataSet);
  }
}

void UnifiedDriverWalkerControlMPITest::reportWalkersPerRank(Communicate* c, MCPopulation& pop)
{
#if !defined(NDEBUG)
  std::vector<int> rank_walker_count(c->size(), 0);
  rank_walker_count[c->rank()] = pop.get_num_local_walkers();
  c->allreduce(rank_walker_count);

  const int current_population = std::accumulate(rank_walker_count.begin(), rank_walker_count.end(), 0);

  if (c->rank() == 0)
  {
    std::cout << "Walkers Per Rank (Total: " << current_population << ")\n";
    for (int i = 0; i < rank_walker_count.size(); ++i)
    {
      std::cout << " " << i << "  " << rank_walker_count[i] << '\n';
    }
  }
#endif
}

void UnifiedDriverWalkerControlMPITest::testNewDistribution(std::vector<int>& minus, std::vector<int>& plus)
{
  int num_ranks = dpools_.comm->size();
  std::vector<int> num_per_rank(num_ranks, 1);
  num_per_rank[0] = num_ranks;
  std::cout << "num per rank:" << NativePrint(num_per_rank) << '\n';
  if (dpools_.comm->rank() == 0)
  {
    auto& walker        = pop_->getWalkerElementsRef(0).walker;
    walker.Multiplicity = num_ranks;
  }
  std::vector<int> fair_offset;
  wc_.determineNewWalkerPopulation(num_per_rank, fair_offset, minus, plus);
}

} // namespace testing

TEST_CASE("WalkerControl::determineNewWalkerPopulation", "[drivers][walker_control]")
{
  std::vector<int> minus;
  std::vector<int> plus;
  testing::UnifiedDriverWalkerControlMPITest test;

  test.testNewDistribution(minus, plus);
  int rank = test.getRank();
  std::cout << "rank:" << rank << " minus: " << NativePrint(minus) << '\n';

  std::cout << "rank:" << rank << " plus: " << NativePrint(plus) << '\n';
  // CHECK(minus.size() == 2);
  // CHECK(plus.size() == 2);
}


void testing::UnifiedDriverWalkerControlMPITest::testPopulationDiff(std::vector<int>& rank_counts_before,
                                                                    std::vector<int>& rank_counts_after)
{
  using MCPWalker = MCPopulation::MCPWalker;

  int rank = dpools_.comm->rank();

  //for (int iw = 0; iw < rank_counts_before
  pop_->get_walkers()[0]->Multiplicity = rank_counts_before[rank];

  std::vector<int> fair_offset;
  std::vector<int> minus;
  std::vector<int> plus;

  wc_.setNumPerRank(rank_counts_before);
  // note this ordering of calls is different from the main code which
  // preserves the ordering to allow the "optimization" of sending multiple copies of the same walker to another rank.
  // I think it should follow this ordering.
  // I only avoids going wrong in the main coded by the fact that
  pop_->copyHighMultiplicityWalkers();
  wc_.swapWalkersSimple(*pop_);
  wc_.killDeadWalkersOnRank(*pop_);
  reportWalkersPerRank(dpools_.comm, *pop_);
  CHECK(pop_->get_num_local_walkers() == rank_counts_after[rank]);
}

void testing::UnifiedDriverWalkerControlMPITest::testWalkerIDs(std::vector<std::vector<int>> walker_ids_after)
{
  int rank = dpools_.comm->rank();
  std::vector<int> walker_ids;
  for (int iw = 0; iw < pop_->get_walkers().size(); ++iw)
    walker_ids.push_back(pop_->get_walkers()[iw]->getWalkerID());
  std::cout << "rank: " << rank << "  walker ids: " << NativePrint(walker_ids) << '\n';

  for (int iw = 0; iw < walker_ids_after[rank].size(); ++iw)
  {
    CHECK(pop_->get_walkers()[iw]->getWalkerID() == walker_ids_after[rank][iw]);
  }
}

TEST_CASE("MPI WalkerControl population swap walkers", "[drivers][walker_control]")
{
  auto test_func = []() {
    outputManager.pause();
    testing::UnifiedDriverWalkerControlMPITest test;
    outputManager.resume();

    SECTION("Simple")
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
      test.testPopulationDiff(count_before, count_after);
      std::vector<long> rank_0_ids(num_ranks);
      std::generate(rank_0_ids.begin(), rank_0_ids.end(),
                    [n = 0, num_ranks]() mutable { return (n++) * num_ranks + 1; });
      std::vector<std::vector<int>> ar_wids;
      for (int ir = 0; ir < num_ranks; ++ir)
      {
        ar_wids.push_back({ir + 1});
        if (ir > 0)
          ar_wids.back().push_back(rank_0_ids[num_ranks - ir]);
      }
      test.testWalkerIDs(ar_wids);
    }
  };
  MPIExceptionWrapper mew;
  mew(test_func);
}

} // namespace qmcplusplus
