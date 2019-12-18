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
    pop_ = std::make_unique<MCPopulation>(num_ranks, dpools_.particle_pool->getParticleSet("e"),
                                         dpools_.wavefunction_pool->getPrimary(),
                                         dpools_.hamiltonian_pool->getPrimary());
    pop_->createWalkers(1);

    wc_.use_nonblocking = true;

    // Set up Cur_pop
    wc_.Cur_pop = dpools_.comm->size();
    for (int i = 0; i < dpools_.comm->size(); i++)
    {
      wc_.NumPerNode[i] = 1;
    }
  }
                                        
  void testMultiplicity()
  {
    using MCPWalker = MCPopulation::MCPWalker;
    // One walker on every node, should be no swapping

    WalkerControlBase::PopulationAdjustment pop_adjust{1,
                                                       convertUPtrToRefVector(pop_->get_walkers()),
                                                       {0},
                                                       RefVector<MCPWalker>{}};
    wc_.swapWalkersSimple(*pop_, pop_adjust);

    REQUIRE(pop_->get_num_local_walkers() == 1);

    // Now bump the multiplicity on rank 0 by 2
    if (dpools_.comm->rank() == 0)
    {
      pop_adjust.copies_to_make[0] += 2;
    }
    wc_.NumPerNode[0] = 3;
    wc_.Cur_pop += 2;
    reportWalkersPerRank(dpools_.comm, *pop_, pop_adjust);
    wc_.swapWalkersSimple(*pop_, pop_adjust);
    reportWalkersPerRank(dpools_.comm, *pop_, pop_adjust);

  }
  
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
  void reportWalkersPerRank(Communicate* c, MCPopulation& pop, WalkerControlBase::PopulationAdjustment& adjust)
  {
    std::vector<int> rank_walker_count(c->size(), 0);
    std::vector<int> adjust_walker_count(c->size(), 0);
    rank_walker_count[c->rank()] = pop.get_num_local_walkers();
    c->allreduce(rank_walker_count);
    int good_walker_count = 0;
    for (auto copies : adjust.copies_to_make)
      good_walker_count += copies + 1;
    adjust_walker_count[c->rank()] = good_walker_count;
    c->allreduce(adjust_walker_count);

    if (c->rank() == 0)
    {
      std::cout << "Walkers Per Rank (Total: " << wc_.Cur_pop << ")\n";
      for (int i = 0; i < rank_walker_count.size(); ++i)
      {
        std::cout << " " << i << "  " << rank_walker_count[i] << " : " << adjust_walker_count[i] << "\n";
      }
    }
  }

  SetupPools dpools_;
  UPtr<MCPopulation> pop_;
  WalkerControlMPI wc_;
};
} // namespace testing

TEST_CASE("MPI WalkerControl multiplicity swap walkers", "[drivers][walker_control]")
{
  testing::UnifiedDriverWalkerControlMPITest test;
  test.testMultiplicity();
}

// TEST_CASE("MPI Walker Unified Driver swap walkers", "[drivers][walker_control]")
// {
//   testing::UnifiedDriverWalkerControlMPITest test;
//   test();
// }

} // namespace qmcplusplus
