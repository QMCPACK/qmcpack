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
#include "Lattice/Uniform3DGridLayout.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "Particle/SymmetricDistanceTableData.h"
#include "QMCApp/ParticleSetPool.h"
#include "QMCHamiltonians/CoulombPBCAA.h"


#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{

TEST_CASE("Coulomb PBC A-A", "[hamiltonian]")
{

  LRCoulombSingleton::CoulombHandler = 0;

  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  Uniform3DGridLayout grid;
  grid.BoxBConds = true; // periodic
  grid.R.diagonal(1.0);
  grid.reset();


  ParticleSet ions;
  ParticleSet elec;

  ions.setName("ion");
  ions.create(1);
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.Lattice.copy(grid);

  SpeciesSet &ion_species =  ions.getSpeciesSet();
  int pIdx = ion_species.addSpecies("H");
  int pChargeIdx = ion_species.addAttribute("charge");
  int pMembersizeIdx = ion_species.addAttribute("membersize");
  ion_species(pChargeIdx, pIdx) = 1;
  ion_species(pMembersizeIdx, pIdx) = 1;
  ions.Lattice.copy(grid);
  ions.createSK();


  CoulombPBCAA caa = CoulombPBCAA(ions, false);

  // Background charge term 
  double consts = caa.evalConsts();
  REQUIRE(consts == Approx(-3.1151210154));

  double val = caa.evaluate(ions);
  //cout << "val = " << val << std::endl;
  REQUIRE(val == Approx(-1.418648723)); // not validated


}

TEST_CASE("Coulomb PBC A-A BCC H", "[hamiltonian]")
{

  LRCoulombSingleton::CoulombHandler = 0;

  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  Uniform3DGridLayout grid;
  grid.BoxBConds = true; // periodic
  grid.R.diagonal(3.77945227);
  grid.reset();


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
  ions.Lattice.copy(grid);

  SpeciesSet &ion_species =  ions.getSpeciesSet();
  int pIdx = ion_species.addSpecies("H");
  int pChargeIdx = ion_species.addAttribute("charge");
  int pMembersizeIdx = ion_species.addAttribute("membersize");
  ion_species(pChargeIdx, pIdx) = 1;
  ion_species(pMembersizeIdx, pIdx) = 2;
  ions.Lattice.copy(grid);
  ions.createSK();


  CoulombPBCAA caa = CoulombPBCAA(ions, false);

  // Background charge term 
  double consts = caa.evalConsts();
  REQUIRE(consts == Approx(-1.675229452)); // not validated

  double val = caa.evaluate(elec);
  //cout << "BCC H val = " << val << std::endl;
  REQUIRE(val == Approx(-0.9628996199)); // not validated


}

TEST_CASE("Coulomb PBC A-A elec", "[hamiltonian]")
{

  LRCoulombSingleton::CoulombHandler = 0;

  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  Uniform3DGridLayout grid;
  grid.BoxBConds = true; // periodic
  grid.R.diagonal(1.0);
  grid.reset();

  ParticleSet elec;

  elec.Lattice.copy(grid);
  elec.setName("elec");
  elec.create(1);
  elec.R[0][0] = 0.0;
  elec.R[0][1] = 0.5;
  elec.R[0][2] = 0.0;

  SpeciesSet &tspecies =  elec.getSpeciesSet();
  int upIdx = tspecies.addSpecies("u");
  int downIdx = tspecies.addSpecies("d");
  int chargeIdx = tspecies.addAttribute("charge");
  int massIdx = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(chargeIdx, downIdx) = -1;
  tspecies(massIdx, upIdx) = 1.0;
  tspecies(massIdx, downIdx) = 1.0;

  elec.createSK();
  elec.update();



  CoulombPBCAA caa = CoulombPBCAA(elec, false);

  // Self-energy correction, no background charge for e-e interaction
  double consts = caa.evalConsts();
  REQUIRE(consts == Approx(-3.064495247));

  double val = caa.evaluate(elec);
  //cout << "val = " << val << std::endl;
  REQUIRE(val == Approx(-1.3680247006)); // not validated


}


}

