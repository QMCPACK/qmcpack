//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "MinimalHamiltonianPool.h"
#include "OhmmsData/Libxml2Doc.h"

namespace qmcplusplus
{

// See src/QMCHamiltonians/tests/test_hamiltonian_factory for parsing tests
static constexpr const char* const hamiltonian_xml = R"(
<hamiltonian name="h0" type="generic" target="e">
  <pairpot type="coulomb" name="ElecElec" source="e" target="e"/>
</hamiltonian>
  )";

static constexpr const char* const hamiltonian_eeei_xml = R"(
<hamiltonian name="h0" type="generic" target="e">
  <pairpot type="coulomb" name="ElecElec" source="e" target="e"/>
  <pairpot type="coulomb" name="ElecIon" source="ion" target="e"/>
</hamiltonian>
  )";

HamiltonianPool MinimalHamiltonianPool::make_hamWithEE(Communicate* comm,
                                                       ParticleSetPool& particle_pool,
                                                       WaveFunctionPool& wavefunction_pool)
{
  HamiltonianPool hpool(particle_pool, wavefunction_pool, comm);
  Libxml2Document doc;
  doc.parseFromString(hamiltonian_xml);

  xmlNodePtr root = doc.getRoot();
  hpool.put(root);

  return hpool;
}

HamiltonianPool MinimalHamiltonianPool::makeHamWithEEEI(Communicate* comm,
                                                        ParticleSetPool& particle_pool,
                                                        WaveFunctionPool& wavefunction_pool)
{
  HamiltonianPool hpool(particle_pool, wavefunction_pool, comm);
  Libxml2Document doc;
  doc.parseFromString(hamiltonian_eeei_xml);

  xmlNodePtr root = doc.getRoot();
  hpool.put(root);

  return hpool;
}

HamiltonianPool MinimalHamiltonianPool::makeHamWithEEEIII(Communicate* comm,
                                                          ParticleSetPool& particle_pool,
                                                          WaveFunctionPool& wavefunction_pool)
{
  HamiltonianPool hpool(particle_pool, wavefunction_pool, comm);
  Libxml2Document doc;
  doc.parseFromString(hamiltonian_eeeiii_xml);

  xmlNodePtr root = doc.getRoot();
  hpool.put(root);

  return hpool;
}

HamiltonianPool MinimalHamiltonianPool::makeHamWithEEEIPS(Communicate* comm,
                                                          ParticleSetPool& particle_pool,
                                                          WaveFunctionPool& wavefunction_pool)
{
  HamiltonianPool hpool(particle_pool, wavefunction_pool, comm);
  Libxml2Document doc;
  doc.parseFromString(hamiltonian_eeeips_xml);

  xmlNodePtr root = doc.getRoot();
  hpool.put(root);

  return hpool;
}


} // namespace qmcplusplus
