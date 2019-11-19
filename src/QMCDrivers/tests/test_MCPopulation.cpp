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


#include "catch.hpp"
#include "Configuration.h"
#include "OhmmsPETE/TinyVector.h"
#include "QMCDrivers/MCPopulation.h"
#include "QMCDrivers/tests/WalkerConsumer.h"
#include "QMCApp/tests/MinimalParticlePool.h"
#include "QMCApp/tests/MinimalWaveFunctionPool.h"
#include "QMCApp/tests/MinimalHamiltonianPool.h"
namespace qmcplusplus
{
TEST_CASE("MCPopulation::createWalkers", "[particle][population]")
{
  using namespace testing;
  Communicate* comm;
  OHMMS::Controller->initialize(0, NULL);
  comm = OHMMS::Controller;

  MinimalParticlePool mpp;
  ParticleSetPool particle_pool = mpp(comm);
  MinimalWaveFunctionPool wfp;
  WaveFunctionPool wavefunction_pool = wfp(comm, &particle_pool);
  wavefunction_pool.setPrimary(wavefunction_pool.getWaveFunction("psi0"));
  MinimalHamiltonianPool mhp;
  HamiltonianPool hamiltonian_pool = mhp(comm, &particle_pool, &wavefunction_pool);

  TrialWaveFunction twf(comm);
  MCPopulation population(1, particle_pool.getParticleSet("e"), &twf, hamiltonian_pool.getPrimary());

  population.createWalkers(8);
  REQUIRE(population.get_walkers().size() == 8);
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

TEST_CASE("MCPopulation::distributeWalkers", "[particle][population]")
{
  using namespace testing;
  Communicate* comm;
  OHMMS::Controller->initialize(0, NULL);
  comm = OHMMS::Controller;

  MinimalParticlePool mpp;
  ParticleSetPool particle_pool = mpp(comm);
  MinimalWaveFunctionPool wfp;
  WaveFunctionPool wavefunction_pool = wfp(comm, &particle_pool);
  wavefunction_pool.setPrimary(wavefunction_pool.getWaveFunction("psi0"));
  MinimalHamiltonianPool mhp;
  HamiltonianPool hamiltonian_pool = mhp(comm, &particle_pool, &wavefunction_pool);

  MCPopulation population(1, particle_pool.getParticleSet("e"), wavefunction_pool.getPrimary(),
                          hamiltonian_pool.getPrimary());

  population.createWalkers(24);
  REQUIRE(population.get_walkers().size() == 24);

  std::vector<std::unique_ptr<WalkerConsumer>> walker_consumers(8);
  std::for_each(walker_consumers.begin(), walker_consumers.end(),
                [](std::unique_ptr<WalkerConsumer>& wc) { wc.reset(new WalkerConsumer()); });
  population.distributeWalkers(walker_consumers.begin(), walker_consumers.end(), 3);
  REQUIRE((*walker_consumers[0]).walkers.size() == 3);

  std::vector<std::unique_ptr<WalkerConsumer>> walker_consumers_incommensurate(5);
  std::for_each(walker_consumers_incommensurate.begin(), walker_consumers_incommensurate.end(),
                [](std::unique_ptr<WalkerConsumer>& wc) { wc.reset(new WalkerConsumer()); });

  population.distributeWalkers(walker_consumers_incommensurate.begin(), walker_consumers_incommensurate.end(), 5);
  REQUIRE((*walker_consumers_incommensurate[0]).walkers.size() == 5);
  REQUIRE((*walker_consumers_incommensurate[4]).walkers.size() == 4);
}

} // namespace qmcplusplus
