//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_MINIMALHAMILTONIANPOOL_H
#define QMCPLUSPLUS_MINIMALHAMILTONIANPOOL_H

#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "QMCHamiltonians/HamiltonianPool.h"
#include "Particle/ParticleSetPool.h"

namespace qmcplusplus
{
class MinimalHamiltonianPool
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

  static constexpr const char* const hamiltonian_eeeiii_xml = R"(
<hamiltonian name="h0" type="generic" target="e">
  <pairpot type="coulomb" name="ElecElec" source="e" target="e"/>
  <pairpot type="coulomb" name="ElecIon" source="ion" target="e"/>
  <pairpot type="coulomb" name="IonIon" source="ion" target="ion"/>
</hamiltonian>
  )";

  static constexpr const char* const hamiltonian_eeeips_xml = R"(
<hamiltonian name="h0" type="generic" target="e">
  <pairpot type="coulomb" name="ElecElec" source="e" target="e"/>
  <pairpot type="coulomb" name="ElecIon" source="ion" target="e"/>
  <pairpot type="pseudo" name="PseudoPot" source="ion0" wavefunction="psi0" format="xml">
     <pseudo elementType="C" href="C.BFD.xml"/>
  </pairpot>
</hamiltonian>
  )";

public:
  static HamiltonianPool make_hamWithEE(Communicate* comm,
                                        ParticleSetPool& particle_pool,
                                        WaveFunctionPool& wavefunction_pool);
  static HamiltonianPool makeHamWithEEEI(Communicate* comm,
                                         ParticleSetPool& particle_pool,
                                         WaveFunctionPool& wavefunction_pool);
  static HamiltonianPool makeHamWithEEEIII(Communicate* comm,
                                           ParticleSetPool& particle_pool,
                                           WaveFunctionPool& wavefunction_pool);
  static HamiltonianPool makeHamWithEEEIPS(Communicate* comm,
                                           ParticleSetPool& particle_pool,
                                           WaveFunctionPool& wavefunction_pool);
};

} // namespace qmcplusplus
#endif
