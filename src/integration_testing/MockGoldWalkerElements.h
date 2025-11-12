//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

/** \file
 *  For many integration tests we need a set of
 *  golden walker elements, in order to reduce the amount of boiler
 *  plate using all three minimal pools to set up tests use this
 *  integration testing helper class
 */

#ifndef QMCPLUSPLUS_MOCK_WALKER_ELEMENTS_FOR_ESTIMATOR_TEST_H
#define QMCPLUSPLUS_MOCK_WALKER_ELEMENTS_FOR_ESTIMATOR_TEST_H

#include <MinimalParticlePool.h>
#include <MinimalWaveFunctionPool.h>
#include <MinimalHamiltonianPool.h>
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

/** make walker elements with primary hamiltonian with Coulombic electron
 *  electron interaction and electron ion interaction
 */
MockGoldWalkerElements makeGoldWalkerElementsWithEE(Communicate*, RuntimeOptions run_time_opt);
/** make walker elements with a primary hamiltonian with Coulombic electron
 *   electron interaction and electron ion interaction and ion ion interaction.
 */
MockGoldWalkerElements makeGoldWalkerElementsWithEEEI(Communicate*, RuntimeOptions run_time_opt);
/** make walker elements with a primary hamiltonian with Coulombic electron
 *   electron interaction and electron ion interaction and a pseudo potential component
 */
MockGoldWalkerElements makeGoldWalkerElementsWithEEEIPS(Communicate*, RuntimeOptions run_time_opt);
} // namespace testing
} // namespace qmcplusplus

#endif
