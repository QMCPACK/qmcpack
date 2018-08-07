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
#include "Particle/DistanceTable.h"
#include "Particle/SymmetricDistanceTableData.h"
#include "QMCApp/ParticleSetPool.h"
#include "QMCHamiltonians/ForceChiesaPBCAA.h"
#include "QMCHamiltonians/ForceCeperley.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"


#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{

// PBC case
TEST_CASE("Chiesa Force Ewald3D", "[hamiltonian]")
{

  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  Uniform3DGridLayout grid;
  grid.BoxBConds = true; // periodic
  grid.R.diagonal(5.0);
  grid.reset();


  ParticleSet ions;
  ParticleSet elec;

  ions.setName("ion");
  ions.create(1);
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;

  elec.setName("elec");
  elec.create(2);
  elec.R[0][0] = 0.0;
  elec.R[0][1] = 1.0;
  elec.R[0][2] = 0.0;
  elec.R[1][0] = 0.4;
  elec.R[1][1] = 0.3;
  elec.R[1][2] = 0.0;

  SpeciesSet &tspecies =  elec.getSpeciesSet();
  int upIdx = tspecies.addSpecies("u");
  int downIdx = tspecies.addSpecies("d");
  //int chargeIdx = tspecies.addAttribute("charge");
  int massIdx = tspecies.addAttribute("mass");
  int eChargeIdx = tspecies.addAttribute("charge");
  tspecies(eChargeIdx, upIdx) = -1.0;
  tspecies(eChargeIdx, downIdx) = -1.0;
  //tspecies(chargeIdx, upIdx) = -1;
  //tspecies(chargeIdx, downIdx) = -1;
  tspecies(massIdx, upIdx) = 1.0;
  tspecies(massIdx, downIdx) = 1.0;

  elec.Lattice.copy(grid);
  elec.createSK();

  SpeciesSet &ion_species =  ions.getSpeciesSet();
  int pIdx = ion_species.addSpecies("H");
  int pChargeIdx = ion_species.addAttribute("charge");
  int pMembersizeIdx = ion_species.addAttribute("membersize");
  ion_species(pChargeIdx, pIdx) = 1;
  ion_species(pMembersizeIdx, pIdx) = 1;
  ions.Lattice.copy(grid);
  ions.createSK();

  ions.resetGroups();

  // The call to resetGroups is needed transfer the SpeciesSet
  // settings to the ParticleSet
  elec.resetGroups();

#ifdef ENABLE_SOA
  elec.addTable(ions,DT_SOA);
  ions.addTable(ions,DT_SOA);
#else
  elec.addTable(ions,DT_AOS);
  ions.addTable(ions,DT_AOS);
#endif

  elec.update();
  ions.update();

  ForceChiesaPBCAA force(ions, elec);
  force.addionion = false;
  force.InitMatrix();

  force.evaluate(elec);
  std::cout << " Force = " << force.forces << std::endl;

  // Unvalidated externally
  REQUIRE(force.forces[0][0] == Approx(3.186559306));
  REQUIRE(force.forces[0][1] == Approx(3.352572459));
  REQUIRE(force.forces[0][2] == Approx(0.0));

}

}

