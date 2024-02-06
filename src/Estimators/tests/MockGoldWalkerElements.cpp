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

#include "MockGoldWalkerElements.h"

namespace qmcplusplus
{
namespace testing
{
MockGoldWalkerElements::MockGoldWalkerElements(Communicate* comm,
                                               RuntimeOptions& runtime_opt,
                                               WaveFunctionPoolFactoryFunc wavefunction_pool_fac_func,
                                               HamPoolFactoryFunc ham_pool_fac_func)
    : particle_pool(MinimalParticlePool::make_diamondC_1x1x1(comm)),
      pset_elec(*(particle_pool.getParticleSet("e"))),
      pset_ions(*(particle_pool.getParticleSet("ion"))),
      wavefunction_pool(wavefunction_pool_fac_func(runtime_opt, comm, particle_pool)),
      hamiltonian_pool(ham_pool_fac_func(comm, particle_pool, wavefunction_pool)),
      twf(*(wavefunction_pool.getPrimary())),
      ham(*(hamiltonian_pool.getPrimary()))
{}

MockGoldWalkerElements makeGoldWalkerElementsWithEE(Communicate* comm, RuntimeOptions runtime_opt)
{
  using namespace std::placeholders;
  MockGoldWalkerElements::WaveFunctionPoolFactoryFunc wfp_diamondC =
      std::bind(&MinimalWaveFunctionPool::make_diamondC_1x1x1, _1, _2, _3);
  MockGoldWalkerElements::HamPoolFactoryFunc hamp_ee = std::bind(MinimalHamiltonianPool::make_hamWithEE, _1, _2, _3);
  return MockGoldWalkerElements(comm, runtime_opt, wfp_diamondC, hamp_ee);
} // namespace qmcplusplus

MockGoldWalkerElements makeGoldWalkerElementsWithEEEI(Communicate* comm, RuntimeOptions runtime_opt)
{
  using namespace std::placeholders;
  MockGoldWalkerElements::WaveFunctionPoolFactoryFunc wfp_diamondC =
      std::bind(&MinimalWaveFunctionPool::make_diamondC_1x1x1, _1, _2, _3);
  MockGoldWalkerElements::HamPoolFactoryFunc hamp_ee = std::bind(MinimalHamiltonianPool::makeHamWithEEEI, _1, _2, _3);
  return MockGoldWalkerElements(comm, runtime_opt, wfp_diamondC, hamp_ee);
}

} // namespace testing
} // namespace qmcplusplus
