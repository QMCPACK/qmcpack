//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Yubo "Paul" Yang, yyang173@illinois.edu, University of Illinois at Urbana-Champaign
//
// File created by: Yubo "Paul" Yang, yyang173@illinois.edu, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "LongRange/EwaldHandler3D.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/StressPBC.h"
#include "Utilities/RuntimeOptions.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
// PBC case
TEST_CASE("Stress BCC H Ewald3D", "[hamiltonian]")
{
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.R.diagonal(3.24957306);
  lattice.LR_dim_cutoff = 40;
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet ions(simulation_cell);
  ParticleSet elec(simulation_cell);

  ions.setName("ion");
  ions.create({2});
  ions.R[0]                     = {0.0, 0.0, 0.0};
  ions.R[1]                     = {1.62478653, 1.62478653, 1.62478653};
  SpeciesSet& ion_species       = ions.getSpeciesSet();
  int pIdx                      = ion_species.addSpecies("H");
  int pChargeIdx                = ion_species.addAttribute("charge");
  ion_species(pChargeIdx, pIdx) = 1;
  ions.resetGroups();
  ions.createSK();
  ions.update();

  elec.setName("elec");
  elec.create({2});
  elec.R[0]                  = {0.4, 0.4, 0.4};
  elec.R[1]                  = {2.02478653, 2.02478653, 2.02478653};
  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  // The call to resetGroups is needed transfer the SpeciesSet
  // settings to the ParticleSet
  elec.resetGroups();
  elec.createSK();

  RuntimeOptions runtime_options;
  TrialWaveFunction psi(runtime_options);

  LRCoulombSingleton::CoulombHandler = std::make_unique<EwaldHandler3D>(ions);
  LRCoulombSingleton::CoulombHandler->initBreakup(ions);
  LRCoulombSingleton::CoulombDerivHandler = std::make_unique<EwaldHandler3D>(ions);
  LRCoulombSingleton::CoulombDerivHandler->initBreakup(ions);

  StressPBC est(ions, elec, psi);

  elec.update();
  est.evaluate(elec);

  // i-i = e-e stress is validated against Quantum Espresso's ewald method
  //  they are also double checked using finite-difference
  CHECK(est.getStressEE()(0, 0) == Approx(-0.01087883));
  CHECK(est.getStressEE()(0, 1) == Approx(0.0));
  CHECK(est.getStressEE()(0, 2) == Approx(0.0));
  CHECK(est.getStressEE()(1, 0) == Approx(0.0));
  CHECK(est.getStressEE()(1, 1) == Approx(-0.01087883));
  CHECK(est.getStressEE()(1, 2) == Approx(0.0));
  CHECK(est.getStressEE()(2, 0) == Approx(0.0));
  CHECK(est.getStressEE()(2, 1) == Approx(0.0));
  CHECK(est.getStressEE()(2, 2) == Approx(-0.01087883));

  //Electron-Ion stress diagonal is internally validated using fd.
  CHECK(est.getStressEI()(0, 0) == Approx(-0.00745376));
  //CHECK(est.getStressEI()(0, 1) == Approx(0.0));
  //CHECK(est.getStressEI()(0, 2) == Approx(0.0));
  //CHECK(est.getStressEI()(1, 0) == Approx(0.0));
  CHECK(est.getStressEI()(1, 1) == Approx(-0.00745376));
  //CHECK(est.getStressEI()(1, 2) == Approx(0.0));
  //CHECK(est.getStressEI()(2, 0) == Approx(0.0));
  //CHECK(est.getStressEI()(2, 1) == Approx(0.0));
  CHECK(est.getStressEI()(2, 2) == Approx(-0.00745376));

  LRCoulombSingleton::CoulombHandler.reset(nullptr);

  LRCoulombSingleton::CoulombDerivHandler.reset(nullptr);
}
} // namespace qmcplusplus
