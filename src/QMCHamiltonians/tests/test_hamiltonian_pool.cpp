//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, mdewing@anl.gov Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"


#include "Configuration.h"
#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "QMCHamiltonians/HamiltonianPool.h"
#include "Particle/ParticleSetPool.h"
#include "QMCWaveFunctions/WaveFunctionPool.h"


#include <stdio.h>
#include <string>
#include <sstream>


namespace qmcplusplus
{
extern std::unique_ptr<ParticleSet> createElectronParticleSet(const SimulationCell& simulation_cell);

TEST_CASE("HamiltonianPool", "[qmcapp]")
{
  Communicate* c;
  c = OHMMS::Controller;

  // See src/QMCHamiltonians/tests/test_hamiltonian_factory for parsing tests
  const char* hamiltonian_xml = "<hamiltonian name=\"h0\" type=\"generic\" target=\"e\"> \
         <pairpot type=\"coulomb\" name=\"ElecElec\" source=\"e\" target=\"e\"/> \
</hamiltonian>";

  Libxml2Document doc;
  bool okay = doc.parseFromString(hamiltonian_xml);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  ParticleSetPool pp(c);
  auto qp = createElectronParticleSet(pp.getSimulationCell());
  pp.addParticleSet(std::move(qp));

  WaveFunctionPool wfp(pp, c);

  WaveFunctionFactory* wf_factory = new WaveFunctionFactory("psi0", *pp.getPool()["e"], pp.getPool(), c);
  wfp.getPool()["psi0"] = wf_factory;
  wfp.setPrimary(wf_factory->getTWF());

  HamiltonianPool hpool(pp, wfp, c);

  REQUIRE(hpool.empty());

  hpool.put(root);

  QMCHamiltonian* h = hpool.getHamiltonian("h0");
  REQUIRE(h != nullptr);

  // Bare kinetic energy is always added
  REQUIRE(h->size() == 2);
}

} // namespace qmcplusplus
