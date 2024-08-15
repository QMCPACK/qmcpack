//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

/** \file
 *  Some ParticleSet functions use the global Random so we need some helper functions to
 *  avoid interminant test state when multiple tests are run from a single test program.
 */

#ifndef QMCPLUSPLUS_MOCK_WALKER_ELEMENTS_FOR_ESTIMATOR_TEST_H
#define QMCPLUSPLUS_MOCK_WALKER_ELEMENTS_FOR_ESTIMATOR_TEST_H

#include "Particle/tests/MinimalParticlePool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
#include "QMCHamiltonians/tests/MinimalHamiltonianPool.h"
#include "Utilities/RuntimeOptions.h"
#include "Message/Communicate.h"
#include <functional>

namespace qmcplusplus
{
namespace testing
{
class MockGoldWalkerElements
{
public:
  using WaveFunctionPoolFactoryFunc =
      std::function<WaveFunctionPool(RuntimeOptions run_time_opt, Communicate* comm, ParticleSetPool& psp)>;
  using HamPoolFactoryFunc =
      std::function<HamiltonianPool(Communicate* comm, ParticleSetPool& psp, WaveFunctionPool& wfp)>;
  MockGoldWalkerElements(Communicate* comm,
                         RuntimeOptions& runtime_opt,
                         WaveFunctionPoolFactoryFunc wfp_func,
                         HamPoolFactoryFunc ham_pool_fac_func);

  ParticleSetPool particle_pool;
  ParticleSet& pset_elec;
  ParticleSet& pset_ions;
  WaveFunctionPool wavefunction_pool;
  HamiltonianPool hamiltonian_pool;
  TrialWaveFunction& twf;
  QMCHamiltonian& ham;
};

MockGoldWalkerElements makeGoldWalkerElementsWithEE(Communicate*, RuntimeOptions run_time_opt);
MockGoldWalkerElements makeGoldWalkerElementsWithEEEI(Communicate*, RuntimeOptions run_time_opt);
} // namespace testing
} // namespace qmcplusplus

#endif
