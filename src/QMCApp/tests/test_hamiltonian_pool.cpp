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
#include "QMCApp/HamiltonianPool.h"
#include "QMCApp/ParticleSetPool.h"
#include "QMCApp/WaveFunctionPool.h"


#include <stdio.h>
#include <string>
#include <sstream>



namespace qmcplusplus
{

ParticleSet *
createElectronParticleSet()
{
  ParticleSet *qp = new ParticleSet;
  qp->setName("e");
  qp->create(2);
  qp->R[0][0] = 1.0;
  qp->R[0][1] = 2.0;
  qp->R[0][2] = 3.0;
  qp->R[1][0] = 0.0;
  qp->R[1][1] = 1.1;
  qp->R[1][2] = 2.2;

  SpeciesSet &tspecies = qp->getSpeciesSet();
  int upIdx = tspecies.addSpecies("u");
  int massIdx = tspecies.addAttribute("mass");
  tspecies(massIdx, upIdx) = 1.0;

  return qp;
}


TEST_CASE("HamiltonianPool", "[qmcapp]")
{

  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  HamiltonianPool hpool(c);

  REQUIRE(hpool.empty());

// See src/QMCHamiltonians/tests/test_hamiltonian_factory for parsing tests
  const char* hamiltonian_xml = \
"<hamiltonian name=\"h0\" type=\"generic\" target=\"e\"> \
         <pairpot type=\"coulomb\" name=\"ElecElec\" source=\"e\" target=\"e\"/> \
</hamiltonian>";

  Libxml2Document doc;
  bool okay = doc.parseFromString(hamiltonian_xml);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  ParticleSetPool pp(c);
  ParticleSet *qp = createElectronParticleSet();
  pp.addParticleSet(qp);

  hpool.setParticleSetPool(&pp);

  WaveFunctionPool wfp(c);
  TrialWaveFunction psi = TrialWaveFunction(c);
  wfp.setParticleSetPool(&pp);
  wfp.setPrimary(&psi);

  WaveFunctionFactory::PtclPoolType ptcl_pool;
  ptcl_pool["e"] = qp;
  WaveFunctionFactory *wf_factory = new WaveFunctionFactory(qp, ptcl_pool, c);
  wf_factory->setPsi(&psi);
  wfp.getPool()["psi0"] = wf_factory;

  hpool.setWaveFunctionPool(&wfp);

  hpool.put(root);

  QMCHamiltonian *h = hpool.getHamiltonian("h0");
  REQUIRE(h != NULL);

  // Bare kinetic energy is always added
  REQUIRE(h->size() == 2);
}

}
