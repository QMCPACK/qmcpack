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
#include "QMCHamiltonians/ForceChiesaPBCAA.h"
#include "QMCHamiltonians/ForceCeperley.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "LongRange/EwaldHandler3D.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
// PBC case
TEST_CASE("Chiesa Force BCC H Ewald3D", "[hamiltonian]")
{
  Communicate* c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> Lattice;
  Lattice.BoxBConds = true; // periodic
  Lattice.R.diagonal(3.77945227);
  Lattice.LR_dim_cutoff = 40;
  Lattice.reset();


  ParticleSet ions;
  ParticleSet elec;

  ions.setName("ion");
  ions.create(2);
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.R[1][0] = 1.6;
  ions.R[1][1] = 1.6;
  ions.R[1][2] = 1.88972614;


  SpeciesSet& ion_species       = ions.getSpeciesSet();
  int pIdx                      = ion_species.addSpecies("H");
  int pChargeIdx                = ion_species.addAttribute("charge");
  ion_species(pChargeIdx, pIdx) = 1;
  ions.Lattice = Lattice;
  ions.createSK();


  elec.Lattice = Lattice;
  elec.setName("elec");
  elec.create(2);
  elec.R[0][0] = 0.5;
  elec.R[0][1] = 0.0;
  elec.R[0][2] = 0.0;
  elec.R[1][0] = 0.0;
  elec.R[1][1] = 0.5;
  elec.R[1][2] = 0.0;


  SpeciesSet& tspecies         = elec.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int downIdx                  = tspecies.addSpecies("d");
  int chargeIdx                = tspecies.addAttribute("charge");
  int massIdx                  = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;
  tspecies(massIdx, upIdx)     = 1.0;
  tspecies(massIdx, downIdx)   = 1.0;

  elec.createSK();

  ParticleSetPool ptcl = ParticleSetPool(c);

  ions.resetGroups();

  // The call to resetGroups is needed transfer the SpeciesSet
  // settings to the ParticleSet
  elec.resetGroups();

  LRCoulombSingleton::CoulombHandler = new EwaldHandler3D(ions);
  LRCoulombSingleton::CoulombHandler->initBreakup(ions);
  LRCoulombSingleton::CoulombDerivHandler = new EwaldHandler3D(ions);
  LRCoulombSingleton::CoulombDerivHandler->initBreakup(ions);

  ForceChiesaPBCAA force(ions, elec);
  force.addionion = false;
  force.InitMatrix();

  elec.update();
  force.evaluate(elec);

  //Ion-Ion forces are validated against Quantum Espresso's ewald method:
  REQUIRE(force.forces_IonIon[0][0] == Approx(-0.0228366));
  REQUIRE(force.forces_IonIon[0][1] == Approx(-0.0228366));
  REQUIRE(force.forces_IonIon[0][2] == Approx(0.0000000));
  REQUIRE(force.forces_IonIon[1][0] == Approx(0.0228366));
  REQUIRE(force.forces_IonIon[1][1] == Approx(0.0228366));
  REQUIRE(force.forces_IonIon[1][2] == Approx(0.0000000));

  //Electron-Ion forces are unvalidated externally:
  REQUIRE(force.forces[0][0] == Approx(3.959178977));
  REQUIRE(force.forces[0][1] == Approx(3.959178977));
  REQUIRE(force.forces[0][2] == Approx(0.000000000));
  REQUIRE(force.forces[1][0] == Approx(-0.078308730));
  REQUIRE(force.forces[1][1] == Approx(-0.078308730));
  REQUIRE(force.forces[1][2] == Approx(0.000000000));


  delete LRCoulombSingleton::CoulombHandler;
  LRCoulombSingleton::CoulombHandler = 0;

  delete LRCoulombSingleton::CoulombDerivHandler;
  LRCoulombSingleton::CoulombDerivHandler = 0;
}

} // namespace qmcplusplus
