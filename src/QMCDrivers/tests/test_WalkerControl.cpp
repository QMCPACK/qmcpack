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

#include <math.h>

#include "test_WalkerControl.h"
#include "Message/Communicate.h"
#include "ParticleSetPool.h"
#include "WaveFunctionPool.h"
#include "HamiltonianPool.h"
#include "QMCDrivers/MCPopulation.h"
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
  std::vector<int> num_per_rank = {3, 1, 1};
  std::vector<int> fair_offset;
  WalkerControl::determineNewWalkerPopulation(num_per_rank, fair_offset, minus, plus);
}


void UnifiedDriverWalkerControlMPITest::testPopulationDiff(std::vector<int>& rank_counts_before,
                                                           std::vector<int>& rank_counts_after)
{
  using MCPWalker = MCPopulation::MCPWalker;

  int rank = dpools_.comm->rank();

  // first get the walker counts to the expected numbers
  // we do this buy updating the multiplicity and then calling
  pop_->get_walkers()[0]->Weight = rank_counts_before[rank];
  reportWalkersPerRank(dpools_.comm, *pop_);

  wc_.branch(1, *pop_, false);

  reportWalkersPerRank(dpools_.comm, *pop_);
  CHECK(pop_->get_num_local_walkers() == rank_counts_after[rank]);
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

TEST_CASE("MPI WalkerControl population swap walkers", "[drivers][walker_control]")
{
  auto test_func = []() {
    testing::UnifiedDriverWalkerControlMPITest test;

    std::vector<int> count_before{3, 1, 1};
    std::vector<int> count_after{1, 2, 2};
    test.testPopulationDiff(count_before, count_after);
  };
  MPIExceptionWrapper mew;
  mew(test_func);
}

TEST_CASE("WalkerControl::nCopiesOverflowLogic", "[population]")
{
  // less straight forward but important situation
  // this only works if the limit is an exactly representable integer in the platform's double
  double too_large_d = static_cast<double>(std::numeric_limits<long>::max());
  // std::cout << std::setprecision(20) << too_large_d << '\n';
  CHECK(!(std::floor(too_large_d) < static_cast<double>(std::numeric_limits<long>::max())));
  double small_enough_d = nextafter(static_cast<double>(std::numeric_limits<long>::max()), 0.0);
  // std::cout << std::setprecision(20) <<small_enough_d << '\n';
  CHECK(std::floor(small_enough_d) < std::ceil(static_cast<double>(std::numeric_limits<long>::max())));
}

} // NAMESPACE qmcplusplus
