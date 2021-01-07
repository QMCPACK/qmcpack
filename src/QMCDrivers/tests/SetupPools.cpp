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
#include "SetupPools.h"

namespace qmcplusplus
{
namespace testing
{
SetupPools::SetupPools()
{
  comm = OHMMS::Controller;

  std::cout << "For purposes of multithreaded testing max threads is forced to 8" << '\n';
  Concurrency::OverrideMaxCapacity<> override(8);
  
  particle_pool = std::make_unique<ParticleSetPool>(mpp(comm));
  wavefunction_pool = std::make_unique<WaveFunctionPool>(wfp(comm, *particle_pool));
  hamiltonian_pool = std::make_unique<HamiltonianPool>(mhp(comm, *particle_pool, *wavefunction_pool));
}

} // namespace testing
} // namespace qmcplusplus
