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

#ifndef QMCPLUSPLUS_SETUP_DMCTEST_H
#define QMCPLUSPLUS_SETUP_DMCTEST_H

#include "QMCApp/tests/MinimalParticlePool.h"
#include "QMCApp/tests/MinimalWaveFunctionPool.h"
#include "QMCApp/tests/MinimalHamiltonianPool.h"
#include "Concurrency/Info.hpp"
#include "Concurrency/UtilityFunctions.hpp"
#include "QMCDrivers/DMC/DMCBatched.h"

namespace qmcplusplus
{
class SetupDMCTest
{
public:
  SetupDMCTest(int nranks=4) : num_ranks(nranks), qmcdrv_input(3)
  {
    using namespace testing;

    OHMMS::Controller->initialize(0, NULL);
    comm = OHMMS::Controller;

    Concurrency::OverrideMaxThreads<> override(8);

    Libxml2Document doc;
    doc.parseFromString(valid_dmc_input_sections[valid_dmc_input_dmc_batch_index]);
    node = doc.getRoot();

    qmcdrv_input.readXML(node);
    dmcdrv_input.readXML(node);

    particle_pool.reset(new ParticleSetPool(mpp(comm)));
    wavefunction_pool.reset(new WaveFunctionPool(wfp(comm, particle_pool.get())));
    wavefunction_pool->setPrimary(wavefunction_pool->getWaveFunction("psi0"));
    hamiltonian_pool.reset(new HamiltonianPool(mhp(comm, particle_pool.get(), wavefunction_pool.get())));

    if (Concurrency::maxThreads<>() < 8)
      num_crowds = Concurrency::maxThreads<>();
  }

  DMCBatched operator()()
  {
    int num_ranks = comm->size();

    population = MCPopulation(num_ranks, particle_pool->getParticleSet("e"),
                              wavefunction_pool->getPrimary(), hamiltonian_pool->getPrimary());

    QMCDriverInput qmc_input_copy(qmcdrv_input);
    DMCDriverInput dmc_input_copy(dmcdrv_input);
    return {std::move(qmc_input_copy), std::move(dmc_input_copy), population, *(wavefunction_pool->getPrimary()),
          *(hamiltonian_pool->getPrimary()), *(wavefunction_pool), comm};
  }

public:
  MinimalParticlePool mpp;
  MinimalWaveFunctionPool wfp;
  MinimalHamiltonianPool mhp;

  MCPopulation population;

  UPtr<ParticleSetPool> particle_pool;
  UPtr<WaveFunctionPool> wavefunction_pool;
  UPtr<HamiltonianPool> hamiltonian_pool;

  Libxml2Document doc;
  xmlNodePtr node;

  QMCDriverInput qmcdrv_input;
  DMCDriverInput dmcdrv_input;
  int num_ranks;
  int num_crowds = 8;
  Communicate* comm;
};

} // namespace qmcplusplus

#endif
