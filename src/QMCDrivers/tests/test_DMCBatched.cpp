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
#include "QMCDrivers/DMC/DMCDriverInput.h"
#include "QMCDrivers/DMC/DMCBatched.h"
#include "QMCDrivers/tests/ValidQMCInputSections.h"
#include "QMCApp/tests/MinimalParticlePool.h"
#include "QMCApp/tests/MinimalWaveFunctionPool.h"
#include "QMCApp/tests/MinimalHamiltonianPool.h"
#include "Concurrency/Info.hpp"
#include "Concurrency/UtilityFunctions.hpp"

namespace qmcplusplus
{
TEST_CASE("DMCBatched::calc_default_local_walkers", "[drivers]")
{
  using namespace testing;
  Concurrency::OverrideMaxThreads<> override(8);
  Communicate* comm;
  OHMMS::Controller->initialize(0, NULL);
  comm = OHMMS::Controller;

  Libxml2Document doc;
  doc.parseFromString(valid_dmc_input_sections[valid_dmc_input_dmc_batch_index]);
  xmlNodePtr node = doc.getRoot();
  QMCDriverInput qmcdriver_input(3);
  qmcdriver_input.readXML(node);

  MinimalParticlePool mpp;
  ParticleSetPool particle_pool = mpp(comm);
  MinimalWaveFunctionPool wfp(comm);
  WaveFunctionPool wavefunction_pool = wfp(&particle_pool);
  wavefunction_pool.setPrimary(wavefunction_pool.getWaveFunction("psi0"));
  MinimalHamiltonianPool mhp(comm);
  HamiltonianPool hamiltonian_pool = mhp(&particle_pool, &wavefunction_pool);

  int num_ranks  = 4;
  int num_crowds = 8;
  if (Concurrency::maxThreads<>() < 8)
    num_crowds = Concurrency::maxThreads<>();


  auto testWRTWalkersPerRank = [&](int walkers_per_rank) {
    MCPopulation population(num_ranks, particle_pool.getParticleSet("e"), wavefunction_pool.getPrimary(),
                            hamiltonian_pool.getPrimary());
    QMCDriverInput qmcdriver_copy(qmcdriver_input);
    DMCDriverInput dmcdriver_input{};
    DMCBatched dmc_batched(std::move(qmcdriver_copy), std::move(dmcdriver_input), population,
                           *(wavefunction_pool.getPrimary()), *(hamiltonian_pool.getPrimary()), wavefunction_pool,
                           comm);
    dmc_batched.set_walkers_per_rank(walkers_per_rank, "testing");
    if (num_crowds < 8)
      dmc_batched.set_num_crowds(Concurrency::maxThreads(), "Insufficient threads available to match test input");
    DMCBatched::IndexType local_walkers       = dmc_batched.calc_default_local_walkers(walkers_per_rank);
    QMCDriverNew::IndexType walkers_per_crowd = dmc_batched.get_walkers_per_crowd();

    if (walkers_per_rank < num_crowds)
    {
      CHECK(walkers_per_crowd == 1);
      CHECK(local_walkers == num_crowds);
      CHECK(population.get_num_local_walkers() == num_crowds);
      CHECK(population.get_num_global_walkers() == num_crowds * num_ranks);
    }
    else if (walkers_per_rank % num_crowds)
    {
      CHECK(walkers_per_crowd == walkers_per_rank / num_crowds + 1);
      CHECK(local_walkers == walkers_per_crowd * num_crowds);
      CHECK(population.get_num_local_walkers() == walkers_per_crowd * num_crowds);
      CHECK(population.get_num_global_walkers() == walkers_per_crowd * num_crowds * num_ranks);
    }
    else
    {
      CHECK(local_walkers == walkers_per_rank);
      CHECK(population.get_num_local_walkers() == walkers_per_rank);
      CHECK(population.get_num_global_walkers() == walkers_per_rank * num_ranks);
    }
  };
  testWRTWalkersPerRank(7);
  testWRTWalkersPerRank(31);
  testWRTWalkersPerRank(32);
  testWRTWalkersPerRank(33);
}

} // namespace qmcplusplus
