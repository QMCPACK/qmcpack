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

#include <functional>

#include "Message/catch_mpi_main.hpp"

//#include <catch.hpp>

#include "Message/Communicate.h"

#include "QMCDrivers/MCPopulation.h"
#include "QMCDrivers/QMCDriverInput.h"
#include "QMCDrivers/tests/test_WalkerControlMPI.h"
#include "Utilities/OutputManager.h"

//#include "Concurrency/Info.hpp"
//#include "Concurrency/UtilityFunctions.hpp"

namespace qmcplusplus
{
namespace testing
{
UnifiedDriverWalkerControlMPITest::UnifiedDriverWalkerControlMPITest() : wc_(dpools_.comm)
{
  using namespace testing;

  int num_ranks = dpools_.comm->size();
  pop_          = std::make_unique<MCPopulation>(num_ranks, dpools_.particle_pool->getParticleSet("e"),
                                        dpools_.wavefunction_pool->getPrimary(), dpools_.hamiltonian_pool->getPrimary(),
                                        dpools_.comm->rank());

  pop_->createWalkers(1);

  wc_.use_nonblocking = true;


  // Set up Cur_pop
  wc_.Cur_pop = dpools_.comm->size();
  for (int i = 0; i < dpools_.comm->size(); i++)
  {
    wc_.NumPerNode[i] = 1;
  }
}

void UnifiedDriverWalkerControlMPITest::testMultiplicity(std::vector<int>& rank_counts_expanded,
                                                         std::vector<int>& rank_counts_after)
{
  using MCPWalker = MCPopulation::MCPWalker;

  int rank = dpools_.comm->rank();

  // Currently some/all this duplicate state is necessary to have a a successful swap.
  // \todo remove duplicate state effecting population control

  pop_->get_walkers()[0]->Multiplicity = rank_counts_expanded[rank];
  int future_pop                       = std::accumulate(rank_counts_expanded.begin(), rank_counts_expanded.end(), 0);
  WalkerControlBase::PopulationAdjustment pop_adjust{future_pop,
                                                     convertUPtrToRefVector(pop_->get_walkers()),
                                                     {rank_counts_expanded[rank] - 1},
                                                     RefVector<MCPWalker>{}};

  reportWalkersPerRank(dpools_.comm, *pop_);
  wc_.swapWalkersSimple(*pop_, pop_adjust, rank_counts_expanded);
  reportWalkersPerRank(dpools_.comm, *pop_);
  CHECK(pop_->get_num_local_walkers() == rank_counts_after[rank]);
}

void UnifiedDriverWalkerControlMPITest::testPopulationDiff(std::vector<int>& rank_counts_before,
                                                           std::vector<int>& rank_counts_after)
{
  using MCPWalker = MCPopulation::MCPWalker;

  int rank = dpools_.comm->rank();

  pop_->get_walkers()[0]->Multiplicity = rank_counts_before[rank];

  WalkerControlBase::PopulationAdjustment pop_adjust{rank_counts_before[rank],
                                                     convertUPtrToRefVector(pop_->get_walkers()),
                                                     {rank_counts_before[rank] - 1},
                                                     RefVector<MCPWalker>{}};

  wc_.adjustPopulation(*pop_, pop_adjust);
  WalkerControlBase::onRankSpawnKill(*pop_, std::move(pop_adjust));


  wc_.Cur_pop = std::accumulate(rank_counts_before.begin(), rank_counts_before.end(), 0);

  auto proper_number_copies = [](int size) -> std::vector<int> { return std::vector<int>(size, 0); };
  WalkerControlBase::PopulationAdjustment pop_adjust2{rank_counts_before[rank],
                                                      convertUPtrToRefVector(pop_->get_walkers()),
                                                      proper_number_copies(pop_->get_num_local_walkers()),
                                                      RefVector<MCPWalker>{}};

  auto num_per_node = WalkerControlBase::syncFutureWalkersPerRank(dpools_.comm, pop_->get_num_local_walkers());

  reportWalkersPerRank(dpools_.comm, *pop_);
  wc_.swapWalkersSimple(*pop_, pop_adjust2, num_per_node);
  reportWalkersPerRank(dpools_.comm, *pop_);
  CHECK(pop_->get_num_local_walkers() == rank_counts_after[rank]);
}

void UnifiedDriverWalkerControlMPITest::reportWalkersPerRank(Communicate* c, MCPopulation& pop)
{
#if !defined(NDEBUG)
  std::vector<int> rank_walker_count(c->size(), 0);
  rank_walker_count[c->rank()] = pop.get_num_local_walkers();
  c->allreduce(rank_walker_count);

  if (c->rank() == 0)
  {
    std::cout << "Walkers Per Rank (Total: " << wc_.Cur_pop << ")\n";
    for (int i = 0; i < rank_walker_count.size(); ++i)
    {
      std::cout << " " << i << "  " << rank_walker_count[i] << '\n';
    }
  }
#endif
}

} // namespace testing

TEST_CASE("WalkerControlMPI::determineNewWalkerPopulation", "[drivers][walker_control]")
{
  int cur_pop      = 5;
  int num_contexts = 3;

  std::vector<int> num_per_node = {3, 1, 1};

  for (int i = 0; i < num_contexts; ++i)
  {
    std::vector<int> fair_offset;
    std::vector<int> minus;
    std::vector<int> plus;
    WalkerControlMPI::determineNewWalkerPopulation(cur_pop, num_contexts, i, num_per_node, fair_offset, minus, plus);

    CHECK(minus.size() == 2);
    CHECK(plus.size() == 2);
  }
}

/** Here we manipulate just the Multiplicity of a set of 1 walkers per rank
 */
TEST_CASE("MPI WalkerControl multiplicity swap walkers", "[drivers][walker_control]")
{
  outputManager.pause();
  testing::UnifiedDriverWalkerControlMPITest test;
  outputManager.resume();
  SECTION("Simple")
  {
    std::vector<int> count_before{1, 1, 1};
    std::vector<int> count_after{1, 1, 1};

    // One walker on every node, should be no swapping
    test.testMultiplicity(count_before, count_after);
  }

  SECTION("LoadBalance")
  {
    std::vector<int> count_before{3, 1, 1};
    std::vector<int> count_after{1, 2, 2};

    test.testMultiplicity(count_before, count_after);
  }
}

TEST_CASE("MPI WalkerControl population swap walkers", "[drivers][walker_control]")
{
  outputManager.pause();
  testing::UnifiedDriverWalkerControlMPITest test;
  outputManager.resume();

  SECTION("Simple")
  {
    std::vector<int> count_before{1, 1, 1};
    std::vector<int> count_after{1, 1, 1};
    // One walker on every node, should be no swapping
    test.testPopulationDiff(count_before, count_after);
  }

  SECTION("LoadBalance")
  {
    std::vector<int> count_before{3, 1, 1};
    std::vector<int> count_after{1, 2, 2};
    test.testPopulationDiff(count_before, count_after);
  }
}

} // namespace qmcplusplus
