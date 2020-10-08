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
#include "QMCHamiltonians/HamiltonianPool.h"
#include "Particle/ParticleSetPool.h"

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
  MinimalHamiltonianPool() : comm_(nullptr) {}
  HamiltonianPool operator()(Communicate* comm, ParticleSetPool& particle_pool, WaveFunctionPool& wavefunction_pool)
  {
    comm_ = comm;
    HamiltonianPool hpool(particle_pool, wavefunction_pool, comm_);
    Libxml2Document doc;
    doc.parseFromString(hamiltonian_xml);

    xmlNodePtr root = doc.getRoot();
    hpool.put(root);

    return hpool;
  }

private:
  Communicate* comm_;
};

} // namespace qmcplusplus
#endif
