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
#include "QMCApp/tests/MinimalParticlePool.h"
#include "QMCApp/tests/MinimalWaveFunctionPool.h"
#include "QMCApp/tests/MinimalHamiltonianPool.h"
namespace qmcplusplus
{
TEST_CASE("QMCHamiltonian::flex_evaluate", "[hamiltonian]")
{
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

  std::vector<QMCHamiltonian> hamiltonians;
  hamiltonians.emplace_back(*(hamiltonian_pool.getPrimary()));
  hamiltonians.emplace_back(*(hamiltonian_pool.getPrimary()));

  std::vector<ParticleSet> elecs;
  elecs.emplace_back(*(particle_pool.getParticleSet("e")));
  elecs.emplace_back(*(particle_pool.getParticleSet("e")));

  // TODO: finish initializing the elecs.
  //std::vector<QMCHamiltonian::RealType> local_energies(QMCHamiltonian::flex_evaluate(makeRefVector<QMCHamiltonian>(hamiltonians), makeRefVector<ParticleSet>(elecs)));

  //TODO: Would be nice to check some values but I think the system needs a little more setup
}

} // namespace qmcplusplus
