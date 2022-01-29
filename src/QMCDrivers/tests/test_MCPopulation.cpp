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
#include "Configuration.h"
#include "OhmmsPETE/TinyVector.h"
#include "QMCDrivers/MCPopulation.h"
#include "QMCDrivers/tests/WalkerConsumer.h"
#include "Particle/tests/MinimalParticlePool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
#include "QMCHamiltonians/tests/MinimalHamiltonianPool.h"
namespace qmcplusplus
{
TEST_CASE("MCPopulation::createWalkers", "[particle][population]")
{
  using namespace testing;
  Communicate* comm;
  comm = OHMMS::Controller;

  auto particle_pool     = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool = MinimalWaveFunctionPool::make_diamondC_1x1x1(comm, particle_pool);
  wavefunction_pool.setPrimary(wavefunction_pool.getWaveFunction("psi0"));
  auto wf_factory       = wavefunction_pool.getWaveFunctionFactory("wavefunction");
  auto hamiltonian_pool = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);
  TrialWaveFunction twf;
  WalkerConfigurations walker_confs;

  MCPopulation population(1, comm->rank(), walker_confs, particle_pool.getParticleSet("e"), &twf, wf_factory,
                          hamiltonian_pool.getPrimary());

  population.createWalkers(8, 2.0);
  CHECK(population.get_walkers().size() == 8);
  CHECK(population.get_dead_walkers().size() == 8);
  CHECK(population.get_num_local_walkers() == 8);
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
  Communicate* comm;
  comm = OHMMS::Controller;

  auto particle_pool     = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool = MinimalWaveFunctionPool::make_diamondC_1x1x1(comm, particle_pool);
  wavefunction_pool.setPrimary(wavefunction_pool.getWaveFunction("psi0"));
  auto hamiltonian_pool = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);
  auto wf_factory       = wavefunction_pool.getWaveFunctionFactory("wavefunction");
  WalkerConfigurations walker_confs;
  MCPopulation population(1, comm->rank(), walker_confs, particle_pool.getParticleSet("e"),
                          wavefunction_pool.getPrimary(), wf_factory, hamiltonian_pool.getPrimary());

  population.createWalkers(8);
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
