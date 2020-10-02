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

#include <iostream>
#include <catch.hpp>
#include "Concurrency/UtilityFunctions.hpp"
#include "QMCDrivers/tests/SetupPools.h"

namespace qmcplusplus
{
namespace testing
{
SetupPools::SetupPools()
{
  comm = OHMMS::Controller;

  std::cout << "For purposes of multithreaded testing max threads is forced to 8" << '\n';
  Concurrency::OverrideMaxCapacity<> override(8);
  
  particle_pool.reset(new ParticleSetPool(mpp(comm)));
  wavefunction_pool.reset(new WaveFunctionPool(wfp(comm, particle_pool.get())));
  wavefunction_pool->setPrimary(wavefunction_pool->getWaveFunction("psi0"));
  hamiltonian_pool.reset(new HamiltonianPool(mhp(comm, particle_pool.get(), wavefunction_pool.get())));
}

} // namespace testing
} // namespace qmcplusplus
