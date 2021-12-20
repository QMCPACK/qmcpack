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


#include <functional>
#include "catch.hpp"

#include "test_WalkerControl.h"
#include "Message/Communicate.h"
#include "QMCDrivers/MCPopulation.h"
#include "QMCDrivers/QMCDriverInput.h"
#include "Utilities/MPIExceptionWrapper.hpp"
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
  if (num_ranks != 3)
    throw std::runtime_error("Bad Rank Count, WalkerControlMPI tests can only be run with 3 MPI ranks.");
  pop_ = std::make_unique<MCPopulation>(num_ranks, dpools_.comm->rank(), walker_confs,
                                        dpools_.particle_pool->getParticleSet("e"),
                                        dpools_.wavefunction_pool->getPrimary(),
                                        dpools_.wavefunction_pool->getWaveFunctionFactory("wavefunction"),
                                        dpools_.hamiltonian_pool->getPrimary());

  pop_->createWalkers(1);
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
  std::vector<int> num_per_rank = {3, 1, 1};
  std::vector<int> fair_offset;
  WalkerControl::determineNewWalkerPopulation(num_per_rank, fair_offset, minus, plus);
}

} // namespace testing

TEST_CASE("WalkerControl::determineNewWalkerPopulation", "[drivers][walker_control]")
{
  std::vector<int> minus;
  std::vector<int> plus;

  testing::UnifiedDriverWalkerControlMPITest::testNewDistribution(minus, plus);
  CHECK(minus.size() == 2);
  CHECK(plus.size() == 2);
}

/** Here we manipulate just the Multiplicity of a set of 1 walkers per rank
 */
// Fails in debug after PR #2855 run unit tests in debug!
// trips assert in ../src/Particle/Walker.h:514 !

// TEST_CASE("MPI WalkerControl multiplicity swap walkers", "[drivers][walker_control]")
// {
//   auto test_func = []() {
//     outputManager.pause();
//     testing::UnifiedDriverWalkerControlMPITest test;
//     outputManager.resume();
//     test.makeValidWalkers();
//     SECTION("Simple")
//     {
//       std::vector<int> count_before{1, 1, 1};
//       std::vector<int> count_after{1, 1, 1};

//       // One walker on every node, should be no swapping
//       test.testMultiplicity(count_before, count_after);
//     }

//     SECTION("LoadBalance")
//     {
//       std::vector<int> count_before{3, 1, 1};
//       std::vector<int> count_after{1, 2, 2};

//       test.testMultiplicity(count_before, count_after);
//     }
//   };
//   MPIExceptionWrapper mew;
//   mew(test_func);
// }

// Fails in debug after PR #2855 run unit tests in debug!
// trips assert at ../src/Particle/Walker.h:559 !
// TEST_CASE("MPI WalkerControl population swap walkers", "[drivers][walker_control]")
// {
//   auto test_func = []() {
//     outputManager.pause();
//     testing::UnifiedDriverWalkerControlMPITest test;
//     outputManager.resume();

//     SECTION("Simple")
//     {
//       std::vector<int> count_before{1, 1, 1};
//       std::vector<int> count_after{1, 1, 1};
//       // One walker on every node, should be no swapping
//       test.testPopulationDiff(count_before, count_after);
//     }

//     SECTION("LoadBalance")
//     {
//       std::vector<int> count_before{3, 1, 1};
//       std::vector<int> count_after{1, 2, 2};
//       test.testPopulationDiff(count_before, count_after);
//     }
//   };
//   MPIExceptionWrapper mew;
//   mew(test_func);
// }


} // namespace qmcplusplus
