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
struct UnifiedDriverWalkerControlMPITest
{
  void operator()()
  {
    using namespace testing;
    SetupDMCTest dtest;

    WalkerControlMPI wc(dtest.comm);
    // This finishes setup of population as a side effect.
    DMCBatched dmc_batched(std::move(dtest()));

    dtest.population.createWalkers(1);
    
    wc.use_nonblocking = true;

    // Set up Cur_pop
    wc.Cur_pop = dtest.comm->size();
    for (int i = 0; i < dtest.comm->size(); i++)
    {
      wc.NumPerNode[i] = 1;
    }

    using MCPWalker = MCPopulation::MCPWalker;
    // One walker on every node, should be no swapping
    
    WalkerControlBase::PopulationAdjustment pop_adjust{1,
                                                       convertUPtrToRefVector(dtest.population.get_walkers()),
                                                       {1},
                                                       RefVector<MCPWalker>{}};
    wc.swapWalkersSimple(dtest.population, pop_adjust);

    REQUIRE(dtest.population.get_num_local_walkers() == 1);


    std::cout << "Adding 2 walkers on rank 0\n";
    // add two walkers walkers on rank 0
    // update everyones adjust directly since we are not testing pop_adjust here.
    if (dtest.comm->rank() == 0)
    {
      // Use the ID variable to check that the walker content was transmitted
      pop_adjust.good_walkers.push_back(dtest.population.spawnWalker());
      pop_adjust.good_walkers.back().get().ID = dtest.comm->size();
      pop_adjust.copies_to_make.push_back(0);
      pop_adjust.good_walkers.push_back(dtest.population.spawnWalker());
      pop_adjust.good_walkers.back().get().ID = dtest.comm->size() + 1;
      pop_adjust.copies_to_make.push_back(0);
    }
    wc.NumPerNode[0] = 3;
    wc.Cur_pop += 2;
    
    reportWalkersPerRank(dtest.comm, dtest.population);
    
    wc.swapWalkersSimple(dtest.population, pop_adjust);

    //std::cout << " Rank = " << c->rank() << " good size = " << wc.good_w.size() <<
    //          " ID = " << wc.good_w[0]->ID << std::endl;

    reportWalkersPerRank(dtest.comm, dtest.population);

    if (dtest.comm->size() > 1)
    {
      if (dtest.comm->rank() == dtest.comm->size() - 2)
      {
        REQUIRE(dtest.population.get_num_local_walkers() == 2);
        // This check is a bit too restrictive - no guarantee the last walker was the
        //  one transmitted
        // bool okay1 = dtest.population.get_walkers()[1]->ID == dtest.comm->size() ||
        //     dtest.population.get_walkers()[1]->ID == dtest.comm->size() + 1;
        // REQUIRE(okay1);
      }
      else if (dtest.comm->rank() == dtest.comm->size() - 1)
      {
        REQUIRE(dtest.population.get_num_local_walkers() == 2);
        // bool okay2 = dtest.population.get_walkers()[1]->ID == dtest.comm->size() ||
        //     dtest.population.get_walkers()[1]->ID == dtest.comm->size() + 1;
        // REQUIRE(okay2);
      }
      else
      {
        REQUIRE(dtest.population.get_num_local_walkers() == 1);
        //REQUIRE(dtest.population.get_walkers()[0]->ID == dtest.comm->rank());
      }
    }

        // And now the strange case
    // 6 walkers on rank0, 2 on rank1, 2 on rank2
    if (dtest.comm->size() > 2)
    {
      if (dtest.comm->rank() == 0)
      {
        pop_adjust.good_walkers.push_back(dtest.population.spawnWalker());
        pop_adjust.good_walkers.back().get().ID = dtest.comm->size() + 3;
        pop_adjust.copies_to_make.push_back(2);
        pop_adjust.good_walkers.push_back(dtest.population.spawnWalker());
        pop_adjust.good_walkers.back().get().ID = dtest.comm->size() + 2;
        pop_adjust.copies_to_make.push_back(1);
      }
      wc.NumPerNode[0] = 6;
      wc.Cur_pop += 5;

      reportWalkersPerRank(dtest.comm, dtest.population);
      wc.swapWalkersSimple(dtest.population, pop_adjust);
      reportWalkersPerRank(dtest.comm, dtest.population);
      // These are unique walkers
      if (dtest.comm->rank() == dtest.comm->size() - 2)
      {
        CHECK(dtest.population.get_num_local_walkers() == 3);
      }
      else if (dtest.comm->rank() == dtest.comm->size() - 1)
      {
        CHECK(dtest.population.get_num_local_walkers() == 4);
      }
      else
      {
        CHECK(dtest.population.get_num_local_walkers() == 3);
      }
    }
  }
  private:
  void reportWalkersPerRank(Communicate* c, MCPopulation& pop)
  {
    std::vector<int> rank_walker_count(c->size(), 0);
    rank_walker_count[c->rank()] = pop.get_num_local_walkers();
    c->allreduce(rank_walker_count);
    if (c->rank() == 0)
    {
      std::cout << "Walkers Per Rank (Total: " << pop.get_num_global_walkers() << ")\n";
      for (int i = 0; i < rank_walker_count.size(); ++i)
      {
        std::cout << " " << i << "  " << rank_walker_count[i] << "\n";
      }
    }
  }

};
}

TEST_CASE("MPI Walker Unified Driver swap walkers", "[drivers][walker_control]")
{
  testing::UnifiedDriverWalkerControlMPITest test;
  test();
}

} // namespace qmcplusplus
