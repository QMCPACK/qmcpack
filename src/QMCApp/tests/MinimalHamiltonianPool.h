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

#ifndef QMCPLUSPLUS_MINIMALHAMILTONIANPOOL_H
#define QMCPLUSPLUS_MINIMALHAMILTONIANPOOL_H

#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "QMCApp/HamiltonianPool.h"
#include "QMCApp/ParticleSetPool.h"

namespace qmcplusplus
{
class MinimalHamiltonianPool
{
  // See src/QMCHamiltonians/tests/test_hamiltonian_factory for parsing tests
  const char* hamiltonian_xml = R"(
<hamiltonian name="h0" type="generic" target="e"> 
  <pairpot type="coulomb" name="ElecElec" source="e" target="e"/> 
</hamiltonian>
  )";

public:
  MinimalHamiltonianPool(Communicate* c) : comm_(c) {}
  HamiltonianPool operator()(ParticleSetPool* particle_pool, WaveFunctionPool* wavefunction_pool)
  {
    HamiltonianPool hpool(comm_);
    Libxml2Document doc;
    doc.parseFromString(hamiltonian_xml);

    xmlNodePtr root = doc.getRoot();

    hpool.setParticleSetPool(particle_pool);

    hpool.setWaveFunctionPool(wavefunction_pool);
    hpool.put(root);

    return hpool;
  }

private:
  Communicate* comm_;
};

} // namespace qmcplusplus
#endif
