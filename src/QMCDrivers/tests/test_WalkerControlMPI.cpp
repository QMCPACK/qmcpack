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

#include <catch.hpp>

#include "Message/Communicate.h"

#include "QMCDrivers/QMCDriverInput.h"
#include "QMCDrivers/DMC/WalkerControlMPI.h"
#include "QMCDrivers/tests/ValidQMCInputSections.h"
#include "QMCDrivers/tests/SetupDMCTest.h"
#include "QMCApp/tests/MinimalParticlePool.h"
#include "QMCApp/tests/MinimalWaveFunctionPool.h"
#include "QMCApp/tests/MinimalHamiltonianPool.h"

#include "Concurrency/Info.hpp"
#include "Concurrency/UtilityFunctions.hpp"

namespace qmcplusplus
{
namespace testing
{
/** Once there is only one driver type rename
 */
class UnifiedDriverWalkerControlMPITest
{
public:
  UnifiedDriverWalkerControlMPITest() : wc_(dpools_.comm)
  {
    using namespace testing;

    int num_ranks = dpools_.comm->size();
    pop_ =
        std::make_unique<MCPopulation>(num_ranks, dpools_.particle_pool->getParticleSet("e"),
                                       dpools_.wavefunction_pool->getPrimary(), dpools_.hamiltonian_pool->getPrimary());
    pop_->createWalkers(1);

    wc_.use_nonblocking = true;

    // Set up Cur_pop
    wc_.Cur_pop = dpools_.comm->size();
    for (int i = 0; i < dpools_.comm->size(); i++)
    {
      wc_.NumPerNode[i] = 1;
    }
  }

  void testMultiplicity(std::vector<int>& rank_counts_before, std::vector<int>& rank_counts_after)
  {
    using MCPWalker = MCPopulation::MCPWalker;

    int rank = dpools_.comm->rank();

    // Currently some/all this duplicate state is necessary to have a a successful swap.
    // \todo remove duplicate state effecting population control

    pop_->get_walkers()[0]->Multiplicity = rank_counts_before[rank];
    wc_.Cur_pop                          = std::accumulate(rank_counts_before.begin(), rank_counts_before.end(), 0);
    wc_.NumPerNode                       = rank_counts_before;

    WalkerControlBase::PopulationAdjustment pop_adjust{wc_.NumPerNode[rank],
                                                       convertUPtrToRefVector(pop_->get_walkers()),
                                                       {rank_counts_before[rank] - 1},
                                                       RefVector<MCPWalker>{}};

    reportWalkersPerRank(dpools_.comm, *pop_);
    wc_.swapWalkersSimple(*pop_, pop_adjust, rank_counts_before);
    reportWalkersPerRank(dpools_.comm, *pop_);
    CHECK(pop_->get_num_local_walkers() == rank_counts_after[rank]);
  }

  void testPopulationDiff(std::vector<int>& rank_counts_before, std::vector<int>& rank_counts_after)
  {
    using MCPWalker = MCPopulation::MCPWalker;

    int rank = dpools_.comm->rank();

    // Currently some/all this duplicate state is necessary to have a a successful swap.
    // \todo remove duplicate state effecting population control

    pop_->get_walkers()[0]->Multiplicity = rank_counts_before[rank];

    WalkerControlBase::PopulationAdjustment pop_adjust{wc_.NumPerNode[rank],
                                                       convertUPtrToRefVector(pop_->get_walkers()),
                                                       {rank_counts_before[rank] - 1},
                                                       RefVector<MCPWalker>{}};

    wc_.adjustPopulation(*pop_, pop_adjust);

    auto num_per_node = WalkerControlBase::syncFutureWalkersPerRank(dpools_.comm, pop_adjust.num_walkers);

    wc_.Cur_pop = std::accumulate(rank_counts_before.begin(), rank_counts_before.end(), 0);

    reportWalkersPerRank(dpools_.comm, *pop_);
    wc_.swapWalkersSimple(*pop_, pop_adjust, rank_counts_before);
    reportWalkersPerRank(dpools_.comm, *pop_);
    CHECK(pop_->get_num_local_walkers() == rank_counts_after[rank]);
  }

  //   // Now bump the multiplicity on rank 0 by 2
  //   if (dpools_.comm->rank() == 0)
  //   {
  //     pop_adjust.copies_to_make[0] += 2;
  //   }
  //   wc_.NumPerNode[0] = 3;
  //   wc_.Cur_pop += 2;
  //   reportWalkersPerRank(dpools_.comm, *pop_, pop_adjust);
  //   wc_.swapWalkersSimple(*pop_, pop_adjust);
  //   reportWalkersPerRank(dpools_.comm, *pop_, pop_adjust);

  // }

  // void operator()()
  // {

  //   std::cout << "Adding 2 walkers on rank 0\n";
  //   // add two walkers walkers on rank 0
  //   // update everyones adjust directly since we are not testing pop_adjust here.
  //   if (dtest.comm->rank() == 0)
  //   {
  //     // Use the ID variable to check that the walker content was transmitted
  //     pop_adjust.good_walkers.push_back(dtest.population.spawnWalker());
  //     pop_adjust.good_walkers.back().get().ID = dtest.comm->size();
  //     pop_adjust.copies_to_make.push_back(0);
  //     pop_adjust.good_walkers.push_back(dtest.population.spawnWalker());
  //     pop_adjust.good_walkers.back().get().ID = dtest.comm->size() + 1;
  //     pop_adjust.copies_to_make.push_back(0);
  //   }
  //   wc.NumPerNode[0] = 3;
  //   wc.Cur_pop += 2;

  //   reportWalkersPerRank(dtest.comm, dtest.population, pop_adjust);

  //   wc.swapWalkersSimple(dtest.population, pop_adjust);

  //   //std::cout << " Rank = " << c->rank() << " good size = " << wc.good_w.size() <<
  //   //          " ID = " << wc.good_w[0]->ID << std::endl;

  //   reportWalkersPerRank(dtest.comm, dtest.population, pop_adjust);

  //   if (dtest.comm->size() > 1)
  //   {
  //     if (dtest.comm->rank() == dtest.comm->size() - 2)
  //     {
  //       REQUIRE(dtest.population.get_num_local_walkers() == 2);
  //       // This check is a bit too restrictive - no guarantee the last walker was the
  //       //  one transmitted
  //       // bool okay1 = dtest.population.get_walkers()[1]->ID == dtest.comm->size() ||
  //       //     dtest.population.get_walkers()[1]->ID == dtest.comm->size() + 1;
  //       // REQUIRE(okay1);
  //     }
  //     else if (dtest.comm->rank() == dtest.comm->size() - 1)
  //     {
  //       REQUIRE(dtest.population.get_num_local_walkers() == 2);
  //       // bool okay2 = dtest.population.get_walkers()[1]->ID == dtest.comm->size() ||
  //       //     dtest.population.get_walkers()[1]->ID == dtest.comm->size() + 1;
  //       // REQUIRE(okay2);
  //     }
  //     else
  //     {
  //       REQUIRE(dtest.population.get_num_local_walkers() == 1);
  //       //REQUIRE(dtest.population.get_walkers()[0]->ID == dtest.comm->rank());
  //     }
  //   }

  //   pop_adjust = WalkerControlBase::PopulationAdjustment{1,
  //                                                        convertUPtrToRefVector(dtest.population.get_walkers()),
  //                                                        {0},
  //                                                        RefVector<MCPWalker>{}};

  //   // And now the strange case
  //   // 6 walkers on rank0, 2 on rank1, 2 on rank2
  //   if (dtest.comm->size() > 2)
  //   {
  //     if (dtest.comm->rank() == 0)
  //     {
  //       CHECK(pop_adjust.good_walkers.size() == 1);
  //       CHECK(pop_adjust.copies_to_make.size() == 1);
  //       pop_adjust.good_walkers.push_back(dtest.population.spawnWalker());
  //       pop_adjust.good_walkers.back().get().ID = dtest.comm->size() + 3;
  //       pop_adjust.copies_to_make.push_back(2);
  //       pop_adjust.good_walkers.push_back(dtest.population.spawnWalker());
  //       pop_adjust.good_walkers.back().get().ID = dtest.comm->size() + 2;
  //       pop_adjust.copies_to_make.push_back(1);
  //     }
  //     wc.NumPerNode[0] = 6;
  //     wc.Cur_pop += 5;

  //     reportWalkersPerRank(dtest.comm, dtest.population, pop_adjust);
  //     wc.swapWalkersSimple(dtest.population, pop_adjust);
  //     reportWalkersPerRank(dtest.comm, dtest.population, pop_adjust);
  //     // These are unique walkers
  //     if (dtest.comm->rank() == dtest.comm->size() - 2)
  //     {
  //       CHECK(dtest.population.get_num_local_walkers() == 3);
  //     }
  //     else if (dtest.comm->rank() == dtest.comm->size() - 1)
  //     {
  //       CHECK(dtest.population.get_num_local_walkers() == 4);
  //     }
  //     else
  //     {
  //       CHECK(dtest.population.get_num_local_walkers() == 3);
  //     }
  //   }
  // }

private:
  void reportWalkersPerRank(Communicate* c, MCPopulation& pop)
  {
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
  }

  SetupPools dpools_;
  UPtr<MCPopulation> pop_;
  WalkerControlMPI wc_;
};
} // namespace testing

TEST_CASE("WalkerControlMPI::determineNewWalkerPopulation", "[drivers][walker_control]")
{
  int cur_pop      = 5;
  int num_contexts = 3;

  std::vector<int> num_per_node{1, 1, 1};

  for (int i = 0; i < num_contexts; ++i)
  {
    std::vector<int> fair_offset;
    std::vector<int> minus;
    std::vector<int> plus;
    WalkerControlMPI::determineNewWalkerPopulation(cur_pop, num_contexts, i, num_per_node, fair_offset, minus, plus);
    std::cout << fair_offset[0] << ' ' << fair_offset[1] << '\n' << fair_offset[2] << '\n';

    CHECK(minus.size() == plus.size());
  }
}

/** Here we manipulate just the Multiplicity of a set of 1 walkers per rank
 */
TEST_CASE("MPI WalkerControl multiplicity swap walkers", "[drivers][walker_control]")
{
  testing::UnifiedDriverWalkerControlMPITest test;

  SECTION("Simple")
  {
    std::vector<int> simple_count{1, 1, 1};
    std::vector<int> count_after{1, 1, 1};

    // One walker on every node, should be no swapping
    test.testMultiplicity(simple_count, count_after);
  }

  SECTION("LoadBalance")
  {
    std::vector<int> simple_count{3, 1, 1};
    std::vector<int> count_after{1, 2, 2};

    test.testMultiplicity(simple_count, count_after);
  }
}

TEST_CASE("MPI WalkerControl population swap walkers", "[drivers][walker_control]")
{
  testing::UnifiedDriverWalkerControlMPITest test;

  SECTION("Simple")
  {
    std::vector<int> simple_count{1, 1, 1};
    std::vector<int> count_after{1, 1, 1};
    // One walker on every node, should be no swapping
    test.testPopulationDiff(simple_count, count_after);
  }

  SECTION("LoadBalance")
  {
    std::vector<int> simple_count{3, 1, 1};
    std::vector<int> count_after{1, 2, 2};
    test.testPopulationDiff(simple_count, count_after);
  }
}

// TEST_CASE("MPI Walker Unified Driver swap walkers", "[drivers][walker_control]")
// {
//   testing::UnifiedDriverWalkerControlMPITest test;
//   test();
// }

} // namespace qmcplusplus
