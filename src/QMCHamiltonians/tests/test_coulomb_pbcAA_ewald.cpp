//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Lattice/ParticleBConds.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#ifndef ENABLE_SOA
#include "Particle/SymmetricDistanceTableData.h"
#endif
#include "QMCApp/ParticleSetPool.h"
#include "QMCHamiltonians/CoulombPBCAA.h"
#include "LongRange/EwaldHandler3D.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
TEST_CASE("Coulomb PBC A-A Ewald3D", "[hamiltonian]")
{
  Communicate* c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> Lattice;
  Lattice.BoxBConds = true; // periodic
  Lattice.R.diagonal(1.0);
  Lattice.reset();


  ParticleSet ions;
  ParticleSet elec;

  ions.setName("ion");
  ions.create(1);
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.Lattice = Lattice;

  SpeciesSet& ion_species           = ions.getSpeciesSet();
  int pIdx                          = ion_species.addSpecies("H");
  int pChargeIdx                    = ion_species.addAttribute("charge");
  int pMembersizeIdx                = ion_species.addAttribute("membersize");
  ion_species(pChargeIdx, pIdx)     = 1;
  ion_species(pMembersizeIdx, pIdx) = 1;
  ions.Lattice = Lattice;
  ions.createSK();

  LRCoulombSingleton::CoulombHandler = new EwaldHandler3D(ions);
  LRCoulombSingleton::CoulombHandler->initBreakup(ions);


  CoulombPBCAA caa = CoulombPBCAA(ions, false);
  // Background charge term
  double consts = caa.evalConsts();
  REQUIRE(consts == Approx(-3.142553)); // not validated

  double val = caa.evaluate(ions);
  REQUIRE(val == Approx(-1.418927)); // not validated

  delete LRCoulombSingleton::CoulombHandler;
  LRCoulombSingleton::CoulombHandler = 0;
}

TEST_CASE("Coulomb PBC A-A BCC H Ewald3D", "[hamiltonian]")
{
  Communicate* c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> Lattice;
  Lattice.BoxBConds = true; // periodic
  Lattice.R.diagonal(3.77945227);
  Lattice.reset();


  ParticleSet ions;
  ParticleSet elec;

  ions.setName("ion");
  ions.create(2);
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.R[1][0] = 1.88972614;
  ions.R[1][1] = 1.88972614;
  ions.R[1][2] = 1.88972614;
  ions.Lattice = Lattice;

  SpeciesSet& ion_species           = ions.getSpeciesSet();
  int pIdx                          = ion_species.addSpecies("H");
  int pChargeIdx                    = ion_species.addAttribute("charge");
  int pMembersizeIdx                = ion_species.addAttribute("membersize");
  ion_species(pChargeIdx, pIdx)     = 1;
  ion_species(pMembersizeIdx, pIdx) = 2;
  ions.Lattice = Lattice;
  ions.createSK();

  LRCoulombSingleton::CoulombHandler = new EwaldHandler3D(ions);
  LRCoulombSingleton::CoulombHandler->initBreakup(ions);

  CoulombPBCAA caa = CoulombPBCAA(ions, false);

  // Background charge term
  double consts = caa.evalConsts();
  REQUIRE(consts == Approx(-1.690675)); // not validated

  double val = caa.evaluate(elec);
  REQUIRE(val == Approx(-0.963074)); // not validated

  delete LRCoulombSingleton::CoulombHandler;
  LRCoulombSingleton::CoulombHandler = 0;
}

TEST_CASE("Coulomb PBC A-A elec Ewald3D", "[hamiltonian]")
{
  Communicate* c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> Lattice;
  Lattice.BoxBConds = true; // periodic
  Lattice.R.diagonal(1.0);
  Lattice.reset();

  ParticleSet elec;

  elec.Lattice = Lattice;
  elec.setName("elec");
  elec.create(1);
  elec.R[0][0] = 0.0;
  elec.R[0][1] = 0.5;
  elec.R[0][2] = 0.0;

  SpeciesSet& tspecies         = elec.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int downIdx                  = tspecies.addSpecies("d");
  int chargeIdx                = tspecies.addAttribute("charge");
  int massIdx                  = tspecies.addAttribute("mass");
  int pMembersizeIdx           = tspecies.addAttribute("membersize");
  tspecies(pMembersizeIdx, upIdx)   = 1;
  tspecies(pMembersizeIdx, downIdx) = 0;
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;
  tspecies(massIdx, upIdx)     = 1.0;
  tspecies(massIdx, downIdx)   = 1.0;

  elec.createSK();
  elec.update();

  LRCoulombSingleton::CoulombHandler = new EwaldHandler3D(elec);
  LRCoulombSingleton::CoulombHandler->initBreakup(elec);

  CoulombPBCAA caa = CoulombPBCAA(elec, false);

  // Self-energy correction, no background charge for e-e interaction
  double consts = caa.evalConsts();
  REQUIRE(consts == Approx(-3.142553));

  double val = caa.evaluate(elec);
  REQUIRE(val == Approx(-1.418927)); // not validated

  delete LRCoulombSingleton::CoulombHandler;
  LRCoulombSingleton::CoulombHandler = 0;
}


} // namespace qmcplusplus
