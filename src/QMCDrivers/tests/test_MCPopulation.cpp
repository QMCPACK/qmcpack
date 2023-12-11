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


#include "catch.hpp"
#include <vector>
#include <algorithm>

#include "Configuration.h"
#include "OhmmsPETE/TinyVector.h"
#include "QMCDrivers/MCPopulation.h"
#include "QMCDrivers/tests/WalkerConsumer.h"
#include "Particle/tests/MinimalParticlePool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
#include "QMCHamiltonians/tests/MinimalHamiltonianPool.h"
#include "Utilities/RuntimeOptions.h"

namespace qmcplusplus
{
TEST_CASE("MCPopulation::createWalkers", "[particle][population]")
{
  using namespace testing;

  RuntimeOptions runtime_options;
  Communicate* comm = OHMMS::Controller;

  auto particle_pool     = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool = MinimalWaveFunctionPool::make_diamondC_1x1x1(runtime_options, comm, particle_pool);
  auto hamiltonian_pool  = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);
  TrialWaveFunction twf(runtime_options);
  WalkerConfigurations walker_confs;

  MCPopulation population(1, comm->rank(), particle_pool.getParticleSet("e"), &twf, hamiltonian_pool.getPrimary());

  population.createWalkers(8, walker_confs, 2.0);
  CHECK(population.get_walkers().size() == 8);
  CHECK(population.get_dead_walkers().size() == 8);
  CHECK(population.get_num_local_walkers() == 8);
  population.saveWalkerConfigurations(walker_confs);
  CHECK(walker_confs.getActiveWalkers() == 8);

  MCPopulation population2(1, comm->rank(), particle_pool.getParticleSet("e"), &twf, hamiltonian_pool.getPrimary());
  // keep 3 only configurations.
  walker_confs.resize(3, 0);
  CHECK(walker_confs.getActiveWalkers() == 3);
  auto old_R00 = walker_confs[0]->R[0][0];
  // first and second configurations are both copied from the gold particle set.
  // must be identical.
  CHECK(walker_confs[1]->R[0][0] == old_R00);
  // modify the first configuration
  auto new_R00 = walker_confs[1]->R[0][0] = 0.3;
  population2.createWalkers(8, walker_confs, 1.0);
  CHECK(population2.get_walkers()[0]->R[0][0] == old_R00);
  CHECK(population2.get_walkers()[1]->R[0][0] == new_R00);
  CHECK(population2.get_walkers()[2]->R[0][0] == old_R00);
  CHECK(population2.get_walkers()[3]->R[0][0] == old_R00);
  CHECK(population2.get_walkers()[4]->R[0][0] == new_R00);
  CHECK(population2.get_walkers()[5]->R[0][0] == old_R00);
  CHECK(population2.get_walkers()[6]->R[0][0] == old_R00);
  CHECK(population2.get_walkers()[7]->R[0][0] == new_R00);
}


TEST_CASE("MCPopulation::createWalkers_walker_ids", "[particle][population]")
{
  using namespace testing;

  RuntimeOptions runtime_options;
  Communicate* comm = OHMMS::Controller;

  auto particle_pool     = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool = MinimalWaveFunctionPool::make_diamondC_1x1x1(runtime_options, comm, particle_pool);
  auto hamiltonian_pool  = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);
  TrialWaveFunction twf(runtime_options);
  WalkerConfigurations walker_confs;

  std::vector<MCPopulation> pops;

  int num_ranks = 3;
  for (int i = 0; i < num_ranks; ++i)
    pops.emplace_back(num_ranks, i, particle_pool.getParticleSet("e"), &twf, hamiltonian_pool.getPrimary());

  std::vector<long> walker_ids;
  for (int i = 0; i < num_ranks; ++i)
  {
    pops[i].createWalkers(8, walker_confs, 2.0);
    CHECK(pops[i].get_walkers().size() == 8);
    CHECK(pops[i].get_dead_walkers().size() == 8);
    CHECK(pops[i].get_num_local_walkers() == 8);
    auto walker_elems = pops[i].get_walker_elements();
    for (WalkerElementsRef& wer : walker_elems)
    {
      walker_ids.push_back(wer.walker.ID);
    }
  }
  std::sort(walker_ids.begin(), walker_ids.end());
  // Walker IDs cannot collide
  for (int i = 1; i < walker_ids.size(); ++i)
    CHECK(walker_ids[i - 1] != walker_ids[i]);

  int new_walkers = 3;

  for (int i = 0; i < num_ranks; ++i)
    for (int iw = 0; iw < new_walkers; ++iw)
    {
      auto wer = pops[i].spawnWalker();
      walker_ids.push_back(wer.walker.ID);
    }

  std::sort(walker_ids.begin(), walker_ids.end());
  // Walker IDs cannot collide
  for (int i = 1; i < walker_ids.size(); ++i)
    CHECK(walker_ids[i - 1] != walker_ids[i]);
}


// TEST_CASE("MCPopulation::createWalkers first touch", "[particle][population]")
// {
//   MCPopulation population(1, 2);
//   ParticleAttrib<TinyVector<QMCTraits::RealType, 3>> some_pos(2);
//   some_pos[0] = TinyVector<double, 3>(1.0, 0.0, 0.0);
//   some_pos[1] = TinyVector<double, 3>(0.0, 1.0, 0.0);

//   population.createWalkers(8, 2, 16, some_pos);
//   REQUIRE(population.get_walkers().size() == 16);
//   // Someday here we should use hwloc or similar to query that the memory is close to
//   // each thread.
// }

TEST_CASE("MCPopulation::redistributeWalkers", "[particle][population]")
{
  using namespace testing;

  RuntimeOptions runtime_options;
  Communicate* comm = OHMMS::Controller;

  auto particle_pool     = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool = MinimalWaveFunctionPool::make_diamondC_1x1x1(runtime_options, comm, particle_pool);
  auto hamiltonian_pool  = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);
  WalkerConfigurations walker_confs;
  MCPopulation population(1, comm->rank(), particle_pool.getParticleSet("e"), wavefunction_pool.getPrimary(),
                          hamiltonian_pool.getPrimary());

  population.createWalkers(8, walker_confs);
  REQUIRE(population.get_walkers().size() == 8);

  std::vector<std::unique_ptr<WalkerConsumer>> walker_consumers(2);
  std::for_each(walker_consumers.begin(), walker_consumers.end(),
                [](std::unique_ptr<WalkerConsumer>& wc) { wc.reset(new WalkerConsumer()); });
  population.redistributeWalkers(walker_consumers);

  REQUIRE((*walker_consumers[0]).walkers.size() == 4);

  std::vector<std::unique_ptr<WalkerConsumer>> walker_consumers_incommensurate(3);
  std::for_each(walker_consumers_incommensurate.begin(), walker_consumers_incommensurate.end(),
                [](std::unique_ptr<WalkerConsumer>& wc) { wc.reset(new WalkerConsumer()); });

  population.redistributeWalkers(walker_consumers_incommensurate);
  REQUIRE((*walker_consumers_incommensurate[0]).walkers.size() == 3);
  REQUIRE((*walker_consumers_incommensurate[2]).walkers.size() == 2);
}

} // namespace qmcplusplus
