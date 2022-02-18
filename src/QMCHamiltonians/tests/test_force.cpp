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
#include "QMCHamiltonians/ACForce.h"
#include "QMCHamiltonians/ForceChiesaPBCAA.h"
#include "QMCHamiltonians/ForceCeperley.h"
#include "QMCHamiltonians/CoulombPotential.h"
#include "QMCHamiltonians/CoulombPBCAA.h"
#include "QMCHamiltonians/CoulombPBCAB.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
TEST_CASE("Bare Force", "[hamiltonian]")
{
  Communicate* c = OHMMS::Controller;

  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell);
  ParticleSet elec(simulation_cell);

  ions.setName("ion");
  ions.create({1});
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;

  elec.setName("elec");
  elec.create({2});
  elec.R[0][0] = 0.0;
  elec.R[0][1] = 1.0;
  elec.R[0][2] = 0.0;
  elec.R[1][0] = 0.4;
  elec.R[1][1] = 0.3;
  elec.R[1][2] = 0.0;

  SpeciesSet& tspecies = elec.getSpeciesSet();
  int upIdx            = tspecies.addSpecies("u");
  //int chargeIdx = tspecies.addAttribute("charge");
  int massIdx                 = tspecies.addAttribute("mass");
  int eChargeIdx              = tspecies.addAttribute("charge");
  tspecies(eChargeIdx, upIdx) = -1.0;
  tspecies(massIdx, upIdx)    = 1.0;


  // The call to resetGroups is needed transfer the SpeciesSet
  // settings to the ParticleSet
  elec.resetGroups();

  SpeciesSet& ion_species           = ions.getSpeciesSet();
  int pIdx                          = ion_species.addSpecies("H");
  int pChargeIdx                    = ion_species.addAttribute("charge");
  ion_species(pChargeIdx, pIdx)     = 1;

  ions.resetGroups();
  // Must update ions first in SoA so ions.coordinates_ is valid
  ions.update();

  elec.addTable(ions);
  elec.update();

  BareForce force(ions, elec);
  force.addionion = false;

  force.evaluate(elec);

  //std::cout << " Force = " << force.forces << std::endl;
  REQUIRE(force.forces[0][0] == Approx(3.2));
  REQUIRE(force.forces[0][1] == Approx(3.4));
  REQUIRE(force.forces[0][2] == Approx(0.0));
}

void check_force_copy(ForceChiesaPBCAA& force, ForceChiesaPBCAA& force2)
{
  REQUIRE(force2.Rcut == Approx(force.Rcut));
  REQUIRE(force2.m_exp == force.m_exp);
  REQUIRE(force2.N_basis == force.N_basis);
  REQUIRE(force2.addionion == force.addionion);
  REQUIRE(force2.Sinv.size() == force.Sinv.size());
  std::cout << force.Sinv << std::endl;
  std::cout << force2.Sinv << std::endl;
  for (int i = 0; i < force2.Sinv.rows(); i++)
  {
    for (int j = 0; j < force2.Sinv.cols(); j++)
    {
      //std::cout << "Sinv " << i << "  " << j << " " << force2.Sinv(i,j) << " "  << force.Sinv(i,j) << std::endl;
      REQUIRE(force2.Sinv(i, j) == Approx(force.Sinv(i, j)));
    }
  }

  REQUIRE(force2.h.size() == force.h.size());
  for (int i = 0; i < force2.h.size(); i++)
  {
    REQUIRE(force2.h[i] == Approx(force.h[i]));
  }

  REQUIRE(force2.c.size() == force.c.size());
  for (int i = 0; i < force2.h.size(); i++)
  {
    REQUIRE(force2.c[i] == Approx(force.c[i]));
  }

  REQUIRE(force2.getDistanceTableAAID() == force.getDistanceTableAAID());
  REQUIRE(force2.NumSpeciesA == force.NumSpeciesA);
  REQUIRE(force2.NumSpeciesB == force.NumSpeciesB);
  REQUIRE(force2.NptclA == force.NptclA);
  REQUIRE(force2.NptclB == force.NptclB);
  REQUIRE(force2.Zat.size() == force.Zat.size());
  REQUIRE(force2.Qat.size() == force.Qat.size());
  REQUIRE(force2.Zspec.size() == force.Zspec.size());
  REQUIRE(force2.Qspec.size() == force.Qspec.size());

  REQUIRE(force2.forces_IonIon.size() == force.forces_IonIon.size());
}

// PBC case
TEST_CASE("Chiesa Force", "[hamiltonian]")
{
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.R.diagonal(5.0);
  lattice.LR_dim_cutoff = 25;
  lattice.reset();
  LRCoulombSingleton::this_lr_type = LRCoulombSingleton::EWALD;

  const SimulationCell simulation_cell(lattice);
  ParticleSet ions(simulation_cell);
  ParticleSet elec(simulation_cell);

  ions.setName("ion");
  ions.create({2});
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.R[1][0] = 2.0;
  ions.R[1][1] = 0.0;
  ions.R[1][2] = 0.0;

  elec.setName("elec");
  elec.create({2});
  elec.R[0][0] = 0.0;
  elec.R[0][1] = 1.0;
  elec.R[0][2] = 0.0;
  elec.R[1][0] = 0.4;
  elec.R[1][1] = 0.3;
  elec.R[1][2] = 0.0;

  SpeciesSet& tspecies        = elec.getSpeciesSet();
  int upIdx                   = tspecies.addSpecies("u");
  int massIdx                 = tspecies.addAttribute("mass");
  int eChargeIdx              = tspecies.addAttribute("charge");
  tspecies(eChargeIdx, upIdx) = -1.0;
  tspecies(massIdx, upIdx)    = 1.0;

  elec.createSK();

  SpeciesSet& ion_species           = ions.getSpeciesSet();
  int pIdx                          = ion_species.addSpecies("H");
  int pChargeIdx                    = ion_species.addAttribute("charge");
  ion_species(pChargeIdx, pIdx)     = 1;
  ions.createSK();

  ions.resetGroups();

  // The call to resetGroups is needed transfer the SpeciesSet
  // settings to the ParticleSet
  elec.resetGroups();


  ForceChiesaPBCAA force(ions, elec);
  force.addionion = false;
  force.InitMatrix();

  elec.update();
  force.evaluate(elec);
  std::cout << " Force = " << force.forces << std::endl;
  std::cout << " Forces_IonIon = " << force.forces_IonIon << std::endl;

  // Unvalidated externally
  REQUIRE(force.forces[0][0] == Approx(3.186559306));
  REQUIRE(force.forces[0][1] == Approx(3.352572459));
  REQUIRE(force.forces[0][2] == Approx(0.0));
  REQUIRE(force.forces_IonIon[0][0] == Approx(-0.1478626893));
  REQUIRE(force.forces_IonIon[0][1] == Approx(0.0));
  REQUIRE(force.forces_IonIon[0][2] == Approx(0.0));
  REQUIRE(force.forces_IonIon[1][0] == Approx(0.1478626893));
  REQUIRE(force.forces_IonIon[1][1] == Approx(0.0));
  REQUIRE(force.forces_IonIon[1][2] == Approx(0.0));

  // Let's test CoulombPBCAA and CoulombPBCAB forces, too; Unvalidated externally
  CoulombPBCAA ionForce(ions, false, true, false);
  REQUIRE(ionForce.forces[0][0] == Approx(-0.1478626893));
  REQUIRE(ionForce.forces[1][0] == Approx(0.1478626893));

  CoulombPBCAB elecIonForce(ions, elec, true);
  elecIonForce.evaluate(elec); // Not computed upon construction
  std::cout << " CoulombElecIon = " << elecIonForce.forces << std::endl;
  REQUIRE(elecIonForce.forces[0][0] == Approx(3.186558296));
  REQUIRE(elecIonForce.forces[0][1] == Approx(3.352572459));
  REQUIRE(elecIonForce.forces[1][0] == Approx(-0.3950094326));
  REQUIRE(elecIonForce.forces[1][1] == Approx(0.142639218));

  // The following crafty test is supposed to crash if some checks are out of place
  // This imitates an actual simulation, where Nelec ~= Nnuc that would also crash

  // ParticleSet with 3 ions
  ParticleSet ions3(simulation_cell);
  ions3.setName("ion");
  ions3.create({3});
  ions3.R[0]                           = {0, 0, 0};
  ions3.R[1]                           = {1, 1, 1};
  ions3.R[2]                           = {2, 2, 2};
  SpeciesSet& ion3_species             = ions3.getSpeciesSet();
  int p3Idx                            = ion3_species.addSpecies("H");
  int p3ChargeIdx                      = ion3_species.addAttribute("charge");
  ion3_species(p3ChargeIdx, p3Idx)     = 1;
  ions3.createSK();
  ions3.resetGroups();

  // Namely, sending in incompatible force arrays to evaluateWithIonDerivs is not
  // supposed to do harm, IF
  // 1) forces are not enabled
  CoulombPBCAA noIonForce(ions3, false, false, false);
  // 2) The species is active
  CoulombPBCAA noElecForce(elec, true, true, false);

  TrialWaveFunction psi;
  noIonForce.evaluateWithIonDerivs(ions3, ions3, psi, noElecForce.forces, noElecForce.forces);
  noElecForce.evaluateWithIonDerivs(elec, ions3, psi, noIonForce.forces, noIonForce.forces);

  // It seems a bit silly to test the makeClone method
  // but this class does not use the compiler's copy constructor and
  // there was a bug where the addionion member did not get
  // copied.  Would be nice if there were a better way than inspection
  // to ensure all the members are copied/set up/tested.

  std::unique_ptr<OperatorBase> base_force2 = force.makeClone(elec, psi);
  ForceChiesaPBCAA* force2                  = dynamic_cast<ForceChiesaPBCAA*>(base_force2.get());
  REQUIRE(force2 != nullptr);

  check_force_copy(*force2, force);
}

// Open BC case
TEST_CASE("Ceperley Force", "[hamiltonian]")
{
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell);
  ParticleSet elec(simulation_cell);

  ions.setName("ion");
  ions.create({2});
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.R[1][0] = 2.0;
  ions.R[1][1] = 0.0;
  ions.R[1][2] = 0.0;

  elec.setName("elec");
  elec.create({2});
  elec.R[0][0] = 0.0;
  elec.R[0][1] = 1.0;
  elec.R[0][2] = 0.0;
  elec.R[1][0] = 0.2;
  elec.R[1][1] = 0.3;
  elec.R[1][2] = 0.0;

  SpeciesSet& tspecies        = elec.getSpeciesSet();
  int upIdx                   = tspecies.addSpecies("u");
  int massIdx                 = tspecies.addAttribute("mass");
  int eChargeIdx              = tspecies.addAttribute("charge");
  tspecies(eChargeIdx, upIdx) = -1.0;
  tspecies(massIdx, upIdx)    = 1.0;
  //elec.createSK();

  SpeciesSet& ion_species           = ions.getSpeciesSet();
  int pIdx                          = ion_species.addSpecies("H");
  int pChargeIdx                    = ion_species.addAttribute("charge");
  ion_species(pChargeIdx, pIdx)     = 1;
  //ions.createSK();
  ions.resetGroups();

  // The call to resetGroups is needed transfer the SpeciesSet
  // settings to the ParticleSet
  elec.resetGroups();

  ForceCeperley force(ions, elec);
  force.InitMatrix();

  /// From the 'Force.ipynb' Jupyter notebook
  // for m_exp=2, N_basis=4, Rcut=0.4
  double coeff[4] = {4375, -44296.9, 147656, -161133};
  for (int i = 0; i < 4; i++)
  {
    REQUIRE(force.c[i] == Approx(coeff[i]));
  }

  ions.update();
  elec.update();

  force.addionion = true; // is true by default
  force.evaluate(elec);
  std::cout << " Force ionion = " << force.forces_IonIon << std::endl;
  std::cout << " Force = " << force.forces << std::endl;
  REQUIRE(force.forces[0][0] == Approx(8.99061106).epsilon(1e-4));
  REQUIRE(force.forces[0][1] == Approx(14.86091659).epsilon(1e-4));
  REQUIRE(force.forces[0][2] == Approx(0.0));
  REQUIRE(force.forces[1][0] == Approx(-0.2250998297).epsilon(1e-4));
  REQUIRE(force.forces[1][1] == Approx(0.1388117844).epsilon(1e-4));
  REQUIRE(force.forces[1][2] == Approx(0.0));

  force.N_basis = 6;
  force.Rcut    = 0.8;
  force.InitMatrix();
  // for m_exp=2, N_basis=6, Rcut=0.800000
  double coeff2[6] = {3281.25, -33837.9, 135352, -261841, 245476, -89496.4};
  for (int i = 0; i < 6; i++)
  {
    REQUIRE(force.c[i] == Approx(coeff2[i]));
  }
}

// Test construction of Coulomb forces in OBC
TEST_CASE("Ion-ion Force", "[hamiltonian]")
{
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell);
  ParticleSet elec(simulation_cell);

  ions.setName("ions");
  ions.create({3});
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.R[1][0] = 2.0;
  ions.R[1][1] = 0.0;
  ions.R[1][2] = 0.0;
  ions.R[2][0] = 1.0;
  ions.R[2][1] = 1.0;
  ions.R[2][2] = 0.0;

  // Add elec
  elec.setName("elec");
  elec.create({3});
  elec.R[0][0] = 0.0;
  elec.R[0][1] = 1.0;
  elec.R[0][2] = 0.0;
  elec.R[1][0] = 2.0;
  elec.R[1][1] = 1.0;
  elec.R[1][2] = 0.0;
  elec.R[2][0] = 1.0;
  elec.R[2][1] = 0.0;
  elec.R[2][2] = 0.0;

  SpeciesSet& ionSpecies           = ions.getSpeciesSet();
  int HIdx                         = ionSpecies.addSpecies("H");
  int HChargeIdx                   = ionSpecies.addAttribute("charge");
  ionSpecies(HChargeIdx, HIdx)     = 1;
  ions.resetGroups();

  SpeciesSet& elecSpecies            = elec.getSpeciesSet();
  int upIdx                          = elecSpecies.addSpecies("u");
  int massIdx                        = elecSpecies.addAttribute("mass");
  int eChargeIdx                     = elecSpecies.addAttribute("charge");
  elecSpecies(eChargeIdx, upIdx)     = -1.0;
  elecSpecies(massIdx, upIdx)        = 1.0;
  elec.resetGroups();

  CoulombPotential<OperatorBase::Return_t> ionForce(ions, false, true);
  CoulombPotential<OperatorBase::Return_t> elecIonForce(elec, ions, true); // Should be zero
  CoulombPotential<OperatorBase::Return_t> elecForce(elec, true, true);    // Should be zero

  double coeff0[3] = {-0.60355339059, -0.35355339059, 0.0};
  double coeff1[3] = {0.60355339059, -0.35355339059, 0.0};
  double coeff2[3] = {0.00000000000, 0.70710678119, 0.0};
  for (int i = 0; i < 3; i++)
  {
    REQUIRE(ionForce.forces[0][i] == Approx(coeff0[i]));
    REQUIRE(ionForce.forces[1][i] == Approx(coeff1[i]));
    REQUIRE(ionForce.forces[2][i] == Approx(coeff2[i]));
    REQUIRE(elecIonForce.forces[0][i] == Approx(0.0));
    REQUIRE(elecIonForce.forces[1][i] == Approx(0.0));
    REQUIRE(elecIonForce.forces[2][i] == Approx(0.0));
    REQUIRE(elecForce.forces[0][i] == Approx(0.0));
    REQUIRE(elecForce.forces[1][i] == Approx(0.0));
    REQUIRE(elecForce.forces[2][i] == Approx(0.0));
  }
}

TEST_CASE("AC Force", "[hamiltonian]")
{
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell);
  ParticleSet elec(simulation_cell);

  ions.setName("ion");
  ions.create({1});
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;

  elec.setName("elec");
  elec.create({2});
  elec.R[0][0] = 0.0;
  elec.R[0][1] = 1.0;
  elec.R[0][2] = 0.0;
  elec.R[1][0] = 0.4;
  elec.R[1][1] = 0.3;
  elec.R[1][2] = 0.0;

  SpeciesSet& tspecies = elec.getSpeciesSet();
  int upIdx            = tspecies.addSpecies("u");
  //int chargeIdx = tspecies.addAttribute("charge");
  int massIdx                 = tspecies.addAttribute("mass");
  int eChargeIdx              = tspecies.addAttribute("charge");
  tspecies(eChargeIdx, upIdx) = -1.0;
  tspecies(massIdx, upIdx)    = 1.0;


  // The call to resetGroups is needed transfer the SpeciesSet
  // settings to the ParticleSet
  elec.resetGroups();

  SpeciesSet& ion_species           = ions.getSpeciesSet();
  int pIdx                          = ion_species.addSpecies("H");
  int pChargeIdx                    = ion_species.addAttribute("charge");
  ion_species(pChargeIdx, pIdx)     = 1;

  ions.resetGroups();
  // Must update ions first in SoA so ions.coordinates_ is valid
  ions.update();

  elec.addTable(ions);
  elec.update();

  // defaults
  TrialWaveFunction psi;
  QMCHamiltonian qmcHamiltonian;

  ACForce force(ions, elec, psi, qmcHamiltonian);

  const std::string acforceXML = R"(<tmp> 
  <acforce spacewarp="no" swpow="2." delta="1.e-3">  
  </acforce> 
  </tmp> 
  )";

  Libxml2Document doc;
  bool okay = doc.parseFromString(acforceXML);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();
  xmlNodePtr h1   = xmlFirstElementChild(root);

  force.put(h1);
  const auto v = force.evaluate(elec);
  force.resetTargetParticleSet(elec); // does nothing?

  REQUIRE(v == Approx(0));
  REQUIRE(force.get(std::cout) == true);

  force.add2Hamiltonian(elec, psi, qmcHamiltonian);

  auto clone = force.makeClone(elec, psi, qmcHamiltonian);
  REQUIRE(clone);
}

} // namespace qmcplusplus
