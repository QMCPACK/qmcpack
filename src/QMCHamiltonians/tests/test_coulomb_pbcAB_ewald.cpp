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
#include "QMCHamiltonians/CoulombPBCAB.h"
#include "QMCHamiltonians/CoulombPBCAA.h"
#include "LongRange/EwaldHandler3D.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
TEST_CASE("Coulomb PBC A-B Ewald3D", "[hamiltonian]")
{
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.R.diagonal(1.0);
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet ions(simulation_cell);
  ParticleSet elec(simulation_cell);

  ions.setName("ion");
  ions.create({1});
  ions.R[0]                     = {0.0, 0.0, 0.0};
  SpeciesSet& ion_species       = ions.getSpeciesSet();
  int pIdx                      = ion_species.addSpecies("H");
  int pChargeIdx                = ion_species.addAttribute("charge");
  ion_species(pChargeIdx, pIdx) = 1;
  ions.createSK();
  ions.update();


  elec.setName("elec");
  elec.create({1});
  elec.R[0]                  = {0.5, 0.0, 0.0};
  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  elec.resetGroups();
  elec.createSK();
  elec.addTable(ions);
  elec.update();


  LRCoulombSingleton::CoulombHandler = std::make_unique<EwaldHandler3D>(ions);
  LRCoulombSingleton::CoulombHandler->initBreakup(ions);


  CoulombPBCAB cab(ions, elec);

  // Self energy plus Background charge term
  double consts = cab.evalConsts(elec);
  CHECK(consts == Approx(0.0523598776 * 2)); //not validated

  double val_ei = cab.evaluate(elec);
  CHECK(val_ei == Approx(-0.008302 + 0.0523598776 * 2)); //Not validated

  CoulombPBCAA caa_elec(elec, true, false, false);
  CoulombPBCAA caa_ion(ions, false, false, false);
  double val_ee = caa_elec.evaluate(elec);
  double val_ii = caa_ion.evaluate(ions);
  double sum    = val_ee + val_ii + val_ei;

  CHECK(val_ee == Approx(-1.418927));
  CHECK(val_ii == Approx(-1.418927));
  CHECK(sum == Approx(-2.741436)); // Can be validated via Ewald summation elsewhere
                                   // -2.74136517454081

  LRCoulombSingleton::CoulombHandler.reset(nullptr);
}

TEST_CASE("Coulomb PBC A-B BCC H Ewald3D", "[hamiltonian]")
{
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.R.diagonal(3.77945227);
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet ions(simulation_cell);
  ParticleSet elec(simulation_cell);

  ions.setName("ion");
  ions.create({2});
  ions.R[0]                     = {0.0, 0.0, 0.0};
  ions.R[1]                     = {1.88972614, 1.88972614, 1.88972614};
  SpeciesSet& ion_species       = ions.getSpeciesSet();
  int pIdx                      = ion_species.addSpecies("H");
  int pChargeIdx                = ion_species.addAttribute("charge");
  ion_species(pChargeIdx, pIdx) = 1;
  ions.createSK();
  ions.update();


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

  elec.resetGroups();
  elec.createSK();
  elec.addTable(ions);
  elec.update();


  LRCoulombSingleton::CoulombHandler = std::make_unique<EwaldHandler3D>(ions);
  LRCoulombSingleton::CoulombHandler->initBreakup(ions);

  CoulombPBCAB cab(ions, elec);

  // Background charge term
  double consts = cab.evalConsts(elec);
  CHECK(consts == Approx(0.0277076538 * 4)); //not validated


  double val_ei = cab.evaluate(elec);
  CHECK(val_ei == Approx(-2.223413 + 0.0277076538 * 4)); //Not validated


  CoulombPBCAA caa_elec(elec, false, false, false);
  CoulombPBCAA caa_ion(ions, false, false, false);
  double val_ee = caa_elec.evaluate(elec);
  double val_ii = caa_ion.evaluate(ions);
  double sum    = val_ee + val_ii + val_ei;

  CHECK(val_ee == Approx(-0.012808 - 0.0277076538 * 2));
  CHECK(val_ii == Approx(-0.907659 - 0.0277076538 * 2));
  CHECK(sum == Approx(-3.143880)); // Can be validated via Ewald summation elsewhere
                                   // -3.14349127313640

  LRCoulombSingleton::CoulombHandler.reset(nullptr);
}

} // namespace qmcplusplus
