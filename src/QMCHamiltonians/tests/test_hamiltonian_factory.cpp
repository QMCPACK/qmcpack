//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Configuration.h"
#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "QMCHamiltonians/HamiltonianFactory.h"

namespace qmcplusplus
{
ParticleSet* createElectronParticleSet()
{
  ParticleSet* qp = new ParticleSet;
  qp->setName("e");
  qp->create(2);
  qp->R[0][0] = 1.0;
  qp->R[0][1] = 2.0;
  qp->R[0][2] = 3.0;
  qp->R[1][0] = 0.0;
  qp->R[1][1] = 1.1;
  qp->R[1][2] = 2.2;

  SpeciesSet& tspecies     = qp->getSpeciesSet();
  int upIdx                = tspecies.addSpecies("u");
  int massIdx              = tspecies.addAttribute("mass");
  tspecies(massIdx, upIdx) = 1.0;

  return qp;
}

TEST_CASE("HamiltonianFactory", "[hamiltonian]")
{
  Communicate* c;
  c = OHMMS::Controller;

  ParticleSet* qp = createElectronParticleSet();

  ParticleSet ions;
  ions.setName("ion0");
  ions.create(1);


  HamiltonianFactory::PtclPoolType particle_set_map;
  HamiltonianFactory::PsiPoolType psi_map;

  particle_set_map["e"]    = qp;
  particle_set_map["ion0"] = &ions;

  HamiltonianFactory hf("h0", *qp, particle_set_map, psi_map, c);

  WaveFunctionFactory wff("psi0", *qp, particle_set_map, c);
  psi_map["psi0"] = &wff;

  const char* hamiltonian_xml = "<hamiltonian name=\"h0\" type=\"generic\" target=\"e\"> \
         <pairpot type=\"coulomb\" name=\"ElecElec\" source=\"e\" target=\"e\"/> \
         <pairpot type=\"coulomb\" name=\"IonIon\" source=\"ion0\" target=\"ion0\"/> \
         <pairpot type=\"coulomb\" name=\"ElecIon\" source=\"ion0\" target=\"e\"/> \
</hamiltonian>";

  Libxml2Document doc;
  bool okay = doc.parseFromString(hamiltonian_xml);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();
  hf.put(root);


  REQUIRE(hf.getH());
  REQUIRE(hf.getH()->size() == 3);
  REQUIRE(hf.getH()->total_size() == 3);

  REQUIRE(hf.getH()->getOperatorType("ElecElec") == "coulomb");
  REQUIRE(hf.getH()->getOperatorType("ElecIon") == "coulomb");
}

TEST_CASE("HamiltonianFactory pseudopotential", "[hamiltonian]")
{
  Communicate* c;
  c = OHMMS::Controller;

  ParticleSet* qp = createElectronParticleSet();

  ParticleSet ions;
  ions.setName("ion0");
  std::vector<int> agroup({1});
  ions.create(agroup);

  SpeciesSet& tspecies           = ions.getSpeciesSet();
  int idx                        = tspecies.addSpecies("C");
  int chargeIdx                  = tspecies.addAttribute("charge");
  int atomicNumberIdx            = tspecies.addAttribute("atomicnumber");
  tspecies(chargeIdx, idx)       = 4;
  tspecies(atomicNumberIdx, idx) = 6;


  HamiltonianFactory::PtclPoolType particle_set_map;
  HamiltonianFactory::PsiPoolType psi_map;

  particle_set_map["e"]    = qp;
  particle_set_map["ion0"] = &ions;

  HamiltonianFactory hf("h0", *qp, particle_set_map, psi_map, c);

  WaveFunctionFactory wff("psi0", *qp, particle_set_map, c);
  psi_map["psi0"] = &wff;


  const char* hamilonian_xml = "<hamiltonian name=\"h0\" type=\"generic\" target=\"e\"> \
    <pairpot type=\"pseudo\" name=\"PseudoPot\" source=\"ion0\" wavefunction=\"psi0\" format=\"xml\"> \
        <pseudo elementType=\"C\" href=\"C.BFD.xml\"/> \
     </pairpot> \
</hamiltonian>";

  Libxml2Document doc;
  bool okay = doc.parseFromString(hamilonian_xml);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();
  hf.put(root);
}

} // namespace qmcplusplus
