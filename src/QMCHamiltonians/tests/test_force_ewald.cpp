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
#include "Utilities/RuntimeOptions.h"

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
  ions.R[0]                     = {0.0, 0.0, 0.0};
  ions.R[1]                     = {1.6, 1.6, 1.88972614};
  SpeciesSet& ion_species       = ions.getSpeciesSet();
  int pIdx                      = ion_species.addSpecies("H");
  int pChargeIdx                = ion_species.addAttribute("charge");
  ion_species(pChargeIdx, pIdx) = 1;
  ions.createSK();


  elec.setName("elec");
  elec.create({2});
  elec.R[0]                  = {0.5, 0.0, 0.0};
  elec.R[1]                  = {0.0, 0.5, 0.0};
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
  force.setAddIonIon(false);
  force.InitMatrix();

  elec.update();
  force.evaluate(elec);

  //Ion-Ion forces are validated against Quantum Espresso's ewald method:
  CHECK(force.getForcesIonIon()[0][0] == Approx(-0.0228366));
  CHECK(force.getForcesIonIon()[0][1] == Approx(-0.0228366));
  CHECK(force.getForcesIonIon()[0][2] == Approx(0.0000000));
  CHECK(force.getForcesIonIon()[1][0] == Approx(0.0228366));
  CHECK(force.getForcesIonIon()[1][1] == Approx(0.0228366));
  CHECK(force.getForcesIonIon()[1][2] == Approx(0.0000000));

  //Electron-Ion forces are unvalidated externally:
  CHECK(force.getForces()[0][0] == Approx(3.959178977));
  CHECK(force.getForces()[0][1] == Approx(3.959178977));
  CHECK(force.getForces()[0][2] == Approx(0.000000000));
  CHECK(force.getForces()[1][0] == Approx(-0.078308730));
  CHECK(force.getForces()[1][1] == Approx(-0.078308730));
  CHECK(force.getForces()[1][2] == Approx(0.000000000));

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
  ions.R[0]                     = {0.0, 0.0, 0.0};
  ions.R[1]                     = {1.6, 1.6, 1.88972614};
  SpeciesSet& ion_species       = ions.getSpeciesSet();
  int pIdx                      = ion_species.addSpecies("H");
  int pChargeIdx                = ion_species.addAttribute("charge");
  ion_species(pChargeIdx, pIdx) = 1;
  ions.createSK();


  elec.setName("elec");
  elec.create({2});
  elec.R[0]                  = {0.5, 0.0, 0.0};
  elec.R[1]                  = {0.0, 0.5, 0.0};
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
  force.setAddIonIon(false);
  force.InitMatrix();

  elec.update();
  // test SR part
  force.evaluateSR(elec);
  CHECK(force.getForces()[0][0] == Approx(1.6938118975));
  CHECK(force.getForces()[0][1] == Approx(1.6938118975));
  CHECK(force.getForces()[0][2] == Approx(0.000000000));
  CHECK(force.getForces()[1][0] == Approx(0.000000000));
  CHECK(force.getForces()[1][1] == Approx(0.000000000));
  CHECK(force.getForces()[1][2] == Approx(0.000000000));
  // test LR part
  force.setForces(0);
  force.evaluateLR(elec);
  CHECK(force.getForces()[0][0] == Approx(2.2653670795));
  CHECK(force.getForces()[0][1] == Approx(2.2653670795));
  CHECK(force.getForces()[0][2] == Approx(0.000000000));
  CHECK(force.getForces()[1][0] == Approx(-0.078308730));
  CHECK(force.getForces()[1][1] == Approx(-0.078308730));
  CHECK(force.getForces()[1][2] == Approx(0.000000000));

  // test cloning !!!! makeClone is not testable
  // example call path:
  //  QMCDrivers/CloneManager::makeClones
  //  QMCHamiltonian::makeClone
  //  OperatorBase::add2Hamiltonian -> ForceChiesaPBCAA::makeClone
  RuntimeOptions runtime_options;
  TrialWaveFunction psi(runtime_options);
  std::unique_ptr<ForceChiesaPBCAA> clone(dynamic_cast<ForceChiesaPBCAA*>(force.makeClone(elec, psi).release()));
  clone->evaluate(elec);
  REQUIRE(clone->getAddIonIon() == force.getAddIonIon());
  CHECK(clone->getForcesIonIon()[0][0] == Approx(-0.0228366));
  CHECK(clone->getForcesIonIon()[0][1] == Approx(-0.0228366));
  CHECK(clone->getForcesIonIon()[0][2] == Approx(0.0000000));
  CHECK(clone->getForcesIonIon()[1][0] == Approx(0.0228366));
  CHECK(clone->getForcesIonIon()[1][1] == Approx(0.0228366));
  CHECK(clone->getForcesIonIon()[1][2] == Approx(0.0000000));

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
  ions.R[0]                     = {0.0, 0.0, 0.0};
  ions.R[1]                     = {1.6, 1.6, 1.88972614};
  ions.R[2]                     = {1.4, 0.0, 0.0};
  SpeciesSet& ion_species       = ions.getSpeciesSet();
  int pIdx                      = ion_species.addSpecies("H");
  int pChargeIdx                = ion_species.addAttribute("charge");
  ion_species(pChargeIdx, pIdx) = 1;
  ions.createSK();

  elec.setName("elec");
  elec.create({2});
  elec.R[0]                  = {0.5, 0.0, 0.0};
  elec.R[1]                  = {0.0, 0.5, 0.0};
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
  force.setAddIonIon(false);

  //Ion-Ion forces are validated against Quantum Espresso's ewald method:
  CHECK(force.getForcesIonIon()[0][0] == Approx(-0.37660901));
  CHECK(force.getForcesIonIon()[0][1] == Approx(-0.02283659));
  CHECK(force.getForcesIonIon()[0][2] == Approx(0.0000000));
  CHECK(force.getForcesIonIon()[1][0] == Approx(0.04012282));
  CHECK(force.getForcesIonIon()[1][1] == Approx(0.066670175));
  CHECK(force.getForcesIonIon()[1][2] == Approx(0.0000000));
  CHECK(force.getForcesIonIon()[2][0] == Approx(0.336486185));
  CHECK(force.getForcesIonIon()[2][1] == Approx(-0.04383358));
  CHECK(force.getForcesIonIon()[2][2] == Approx(0.0000000));

  elec.update();
  force.InitMatrix();
  force.evaluate(elec);
  //Electron-Ion forces are unvalidated externally:
  CHECK(force.getForces()[0][0] == Approx(3.959178977));
  CHECK(force.getForces()[0][1] == Approx(3.959178977));
  CHECK(force.getForces()[0][2] == Approx(0.000000000));
  CHECK(force.getForces()[1][0] == Approx(-0.078308730));
  CHECK(force.getForces()[1][1] == Approx(-0.078308730));
  CHECK(force.getForces()[1][2] == Approx(0.000000000));
  CHECK(force.getForces()[2][0] == Approx(-1.4341388802));
  CHECK(force.getForces()[2][1] == Approx(0.1379375923));
  CHECK(force.getForces()[2][2] == Approx(0.000000000));

  LRCoulombSingleton::CoulombHandler.reset(nullptr);
  LRCoulombSingleton::CoulombDerivHandler.reset(nullptr);
}

} // namespace qmcplusplus
