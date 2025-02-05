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

#include "SetupPools.h"
#include <iostream>
#include "Concurrency/UtilityFunctions.hpp"
#include "Particle/tests/MinimalParticlePool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
#include "QMCHamiltonians/tests/MinimalHamiltonianPool.h"
#include "Utilities/ProjectData.h"

namespace qmcplusplus
{
namespace testing
{
SetupPools::SetupPools()
{
  ProjectData test_project("test", ProjectData::DriverVersion::BATCH);
  comm = OHMMS::Controller;

  app_log() << "For purposes of multithreaded testing max threads is forced to 8" << std::endl;
  Concurrency::OverrideMaxCapacity<> override(8);

  particle_pool     = std::make_unique<ParticleSetPool>(MinimalParticlePool::make_diamondC_1x1x1(comm));
  wavefunction_pool = std::make_unique<WaveFunctionPool>(
      MinimalWaveFunctionPool::make_diamondC_1x1x1(test_project.getRuntimeOptions(), comm, *particle_pool));
  hamiltonian_pool = std::make_unique<HamiltonianPool>(
      MinimalHamiltonianPool::make_hamWithEE(comm, *particle_pool, *wavefunction_pool));
}

SetupPools::~SetupPools() = default;

RandomNumberGeneratorPool::RandomNumberGeneratorPool(const size_t num) : rng_pool(num) {}

RandomNumberGeneratorPool::~RandomNumberGeneratorPool() = default;

RefVector<RandomBase<RandomNumberGeneratorPool::FullPrecRealType>> RandomNumberGeneratorPool::getRngRefs()
{
  RefVector<RandomBase<FullPrecRealType>> rng_refs;
  rng_refs.reserve(rng_pool.size());
  for (auto& rng : rng_pool)
    rng_refs.push_back(rng);
  return rng_refs;
}

} // namespace testing
} // namespace qmcplusplus
