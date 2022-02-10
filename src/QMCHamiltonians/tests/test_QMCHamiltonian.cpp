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

#include "type_traits/template_types.hpp"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Particle/tests/MinimalParticlePool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
#include "QMCHamiltonians/tests/MinimalHamiltonianPool.h"
namespace qmcplusplus
{
TEST_CASE("QMCHamiltonian::flex_evaluate", "[hamiltonian]")
{
  Communicate* comm;
  comm = OHMMS::Controller;

  auto particle_pool     = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool = MinimalWaveFunctionPool::make_diamondC_1x1x1(comm, particle_pool);
  wavefunction_pool.setPrimary(wavefunction_pool.getWaveFunction("psi0"));
  auto hamiltonian_pool = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);

  TrialWaveFunction twf;

  std::vector<QMCHamiltonian> hamiltonians;
  //hamiltonians.emplace_back(*(hamiltonian_pool.getPrimary()));
  //hamiltonians.emplace_back(*(hamiltonian_pool.getPrimary()));

  std::vector<ParticleSet> elecs;
  elecs.emplace_back(*(particle_pool.getParticleSet("e")));
  elecs.emplace_back(*(particle_pool.getParticleSet("e")));

  // TODO: finish initializing the elecs.
  //std::vector<QMCHamiltonian::RealType> local_energies(QMCHamiltonian::flex_evaluate(makeRefVector<QMCHamiltonian>(hamiltonians), makeRefVector<ParticleSet>(elecs)));

  //TODO: Would be nice to check some values but I think the system needs a little more setup
}

} // namespace qmcplusplus
