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
#include "Particle/ParticleSet.h"
#include "Particle/ParticleSetPool.h"
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
  Communicate* c = OHMMS::Controller;

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.R.diagonal(3.77945227);
  lattice.LR_dim_cutoff = 40;
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet ions(simulation_cell);
  ParticleSet elec(simulation_cell);

  ions.setName("ion");
  ions.create({2});
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
  ions.createSK();


  elec.setName("elec");
  elec.create({2});
  elec.R[0][0] = 0.5;
  elec.R[0][1] = 0.0;
  elec.R[0][2] = 0.0;
  elec.R[1][0] = 0.0;
  elec.R[1][1] = 0.5;
  elec.R[1][2] = 0.0;


  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  elec.createSK();

  ParticleSetPool ptcl = ParticleSetPool(c);

  ions.resetGroups();

  // The call to resetGroups is needed transfer the SpeciesSet
  // settings to the ParticleSet
  elec.resetGroups();

  LRCoulombSingleton::CoulombHandler = std::make_unique<EwaldHandler3D>(ions);
  LRCoulombSingleton::CoulombHandler->initBreakup(ions);
  LRCoulombSingleton::CoulombDerivHandler = std::make_unique<EwaldHandler3D>(ions);
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

  LRCoulombSingleton::CoulombHandler.reset(nullptr);
  LRCoulombSingleton::CoulombDerivHandler.reset(nullptr);
}

// test SR and LR pieces separately
TEST_CASE("fccz sr lr clone", "[hamiltonian]")
{
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.R.diagonal(3.77945227);
  lattice.LR_dim_cutoff = 40;
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet ions(simulation_cell);
  ParticleSet elec(simulation_cell);

  ions.setName("ion");
  ions.create({2});
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
  ions.createSK();


  elec.setName("elec");
  elec.create({2});
  elec.R[0][0] = 0.5;
  elec.R[0][1] = 0.0;
  elec.R[0][2] = 0.0;
  elec.R[1][0] = 0.0;
  elec.R[1][1] = 0.5;
  elec.R[1][2] = 0.0;

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;
  elec.createSK();

  // The call to resetGroups is needed transfer the SpeciesSet
  // settings to the ParticleSet
  ions.resetGroups();
  elec.resetGroups();

  LRCoulombSingleton::CoulombHandler = std::make_unique<EwaldHandler3D>(ions);
  LRCoulombSingleton::CoulombHandler->initBreakup(ions);
  LRCoulombSingleton::CoulombDerivHandler = std::make_unique<EwaldHandler3D>(ions);
  LRCoulombSingleton::CoulombDerivHandler->initBreakup(ions);

  ForceChiesaPBCAA force(ions, elec);
  force.addionion = false;
  force.InitMatrix();

  elec.update();
  // test SR part
  force.evaluateSR(elec);
  REQUIRE(force.forces[0][0] == Approx(1.6938118975));
  REQUIRE(force.forces[0][1] == Approx(1.6938118975));
  REQUIRE(force.forces[0][2] == Approx(0.000000000));
  REQUIRE(force.forces[1][0] == Approx(0.000000000));
  REQUIRE(force.forces[1][1] == Approx(0.000000000));
  REQUIRE(force.forces[1][2] == Approx(0.000000000));
  // test LR part
  force.forces = 0;
  force.evaluateLR(elec);
  REQUIRE(force.forces[0][0] == Approx(2.2653670795));
  REQUIRE(force.forces[0][1] == Approx(2.2653670795));
  REQUIRE(force.forces[0][2] == Approx(0.000000000));
  REQUIRE(force.forces[1][0] == Approx(-0.078308730));
  REQUIRE(force.forces[1][1] == Approx(-0.078308730));
  REQUIRE(force.forces[1][2] == Approx(0.000000000));

  // test cloning !!!! makeClone is not testable
  // example call path:
  //  QMCDrivers/CloneManager::makeClones
  //  QMCHamiltonian::makeClone
  //  OperatorBase::add2Hamiltonian -> ForceChiesaPBCAA::makeClone
  TrialWaveFunction psi;
  std::unique_ptr<ForceChiesaPBCAA> clone(dynamic_cast<ForceChiesaPBCAA*>(force.makeClone(elec, psi).release()));
  clone->evaluate(elec);
  REQUIRE(clone->addionion == force.addionion);
  REQUIRE(clone->forces_IonIon[0][0] == Approx(-0.0228366));
  REQUIRE(clone->forces_IonIon[0][1] == Approx(-0.0228366));
  REQUIRE(clone->forces_IonIon[0][2] == Approx(0.0000000));
  REQUIRE(clone->forces_IonIon[1][0] == Approx(0.0228366));
  REQUIRE(clone->forces_IonIon[1][1] == Approx(0.0228366));
  REQUIRE(clone->forces_IonIon[1][2] == Approx(0.0000000));

  LRCoulombSingleton::CoulombHandler.reset(nullptr);
  LRCoulombSingleton::CoulombDerivHandler.reset(nullptr);
}

// 3 H atoms randomly distributed in a box
TEST_CASE("fccz h3", "[hamiltonian]")
{
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.R.diagonal(3.77945227);
  lattice.LR_dim_cutoff = 40;
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet ions(simulation_cell);
  ParticleSet elec(simulation_cell);

  ions.setName("ion");
  ions.create({3});
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.R[1][0] = 1.6;
  ions.R[1][1] = 1.6;
  ions.R[1][2] = 1.88972614;
  ions.R[2][0] = 1.4;
  ions.R[2][1] = 0.0;
  ions.R[2][2] = 0.0;

  SpeciesSet& ion_species       = ions.getSpeciesSet();
  int pIdx                      = ion_species.addSpecies("H");
  int pChargeIdx                = ion_species.addAttribute("charge");
  ion_species(pChargeIdx, pIdx) = 1;
  ions.createSK();

  elec.setName("elec");
  elec.create({2});
  elec.R[0][0] = 0.5;
  elec.R[0][1] = 0.0;
  elec.R[0][2] = 0.0;
  elec.R[1][0] = 0.0;
  elec.R[1][1] = 0.5;
  elec.R[1][2] = 0.0;

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;
  elec.createSK();

  // The call to resetGroups is needed transfer the SpeciesSet
  // settings to the ParticleSet
  ions.resetGroups();
  elec.resetGroups();

  LRCoulombSingleton::CoulombHandler = std::make_unique<EwaldHandler3D>(ions);
  LRCoulombSingleton::CoulombHandler->initBreakup(ions);
  LRCoulombSingleton::CoulombDerivHandler = std::make_unique<EwaldHandler3D>(ions);
  LRCoulombSingleton::CoulombDerivHandler->initBreakup(ions);

  ForceChiesaPBCAA force(ions, elec);
  force.addionion = false;

  //Ion-Ion forces are validated against Quantum Espresso's ewald method:
  REQUIRE(force.forces_IonIon[0][0] == Approx(-0.37660901));
  REQUIRE(force.forces_IonIon[0][1] == Approx(-0.02283659));
  REQUIRE(force.forces_IonIon[0][2] == Approx(0.0000000));
  REQUIRE(force.forces_IonIon[1][0] == Approx(0.04012282));
  REQUIRE(force.forces_IonIon[1][1] == Approx(0.066670175));
  REQUIRE(force.forces_IonIon[1][2] == Approx(0.0000000));
  REQUIRE(force.forces_IonIon[2][0] == Approx(0.336486185));
  REQUIRE(force.forces_IonIon[2][1] == Approx(-0.04383358));
  REQUIRE(force.forces_IonIon[2][2] == Approx(0.0000000));

  elec.update();
  force.InitMatrix();
  force.evaluate(elec);
  //Electron-Ion forces are unvalidated externally:
  REQUIRE(force.forces[0][0] == Approx(3.959178977));
  REQUIRE(force.forces[0][1] == Approx(3.959178977));
  REQUIRE(force.forces[0][2] == Approx(0.000000000));
  REQUIRE(force.forces[1][0] == Approx(-0.078308730));
  REQUIRE(force.forces[1][1] == Approx(-0.078308730));
  REQUIRE(force.forces[1][2] == Approx(0.000000000));
  REQUIRE(force.forces[2][0] == Approx(-1.4341388802));
  REQUIRE(force.forces[2][1] == Approx(0.1379375923));
  REQUIRE(force.forces[2][2] == Approx(0.000000000));

  LRCoulombSingleton::CoulombHandler.reset(nullptr);
  LRCoulombSingleton::CoulombDerivHandler.reset(nullptr);
}

} // namespace qmcplusplus
