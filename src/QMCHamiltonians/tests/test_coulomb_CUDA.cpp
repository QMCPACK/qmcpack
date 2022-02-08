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
#include "Particle/MCWalkerConfiguration.h"
#include "QMCHamiltonians/CoulombPBCAB_CUDA.h"
#include "QMCHamiltonians/CoulombPBCAA_CUDA.h"


#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
TEST_CASE("Coulomb PBC A-B CUDA", "[hamiltonian][CUDA]")
{
  LRCoulombSingleton::CoulombHandler = 0;
  LRCoulombSingleton::this_lr_type   = LRCoulombSingleton::ESLER;

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.R.diagonal(1.0);
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet ions(simulation_cell);
  ParticleSet elec(simulation_cell);

  ions.setName("ion");
  ions.create({1});
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;

  SpeciesSet& ion_species       = ions.getSpeciesSet();
  int pIdx                      = ion_species.addSpecies("H");
  int pChargeIdx                = ion_species.addAttribute("charge");
  ion_species(pChargeIdx, pIdx) = 1;
  ions.createSK();
  ions.update();

  elec.setName("elec");
  elec.create({1});
  elec.R[0][0] = 0.0;
  elec.R[0][1] = 0.5;
  elec.R[0][2] = 0.0;

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  elec.createSK();

  elec.addTable(ions);
  elec.update();

  CoulombPBCAB_CUDA cab(ions, elec);

  // Background charge term
  double consts = cab.evalConsts(elec);
  CHECK(consts == Approx(2 * 0.0506238028));

  double val = cab.evaluate(elec);
  //cout << "CUDA val = " << val << std::endl;
  CHECK(val == Approx(-0.005314032183 + 2 * 0.0506238028)); // not validated
}

TEST_CASE("Coulomb PBC AB CUDA BCC H", "[hamiltonian][CUDA]")
{
  LRCoulombSingleton::CoulombHandler = 0;
  LRCoulombSingleton::this_lr_type   = LRCoulombSingleton::ESLER;

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.R.diagonal(3.77945227);
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet ions(simulation_cell);
  MCWalkerConfiguration elec(simulation_cell);

  ions.setName("ion");
  ions.create({2});
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.R[1][0] = 1.88972614;
  ions.R[1][1] = 1.88972614;
  ions.R[1][2] = 1.88972614;

  SpeciesSet& ion_species       = ions.getSpeciesSet();
  int pIdx                      = ion_species.addSpecies("H");
  int pChargeIdx                = ion_species.addAttribute("charge");
  ion_species(pChargeIdx, pIdx) = 1;
  ions.createSK();
  ions.update();

  elec.setName("elec");
  elec.create({2});
  elec.R[0][0] = 0.5;
  elec.R[0][1] = 0.0;
  elec.R[0][2] = 0.0;
  elec.R[1][0] = 0.0;
  elec.R[1][1] = 0.5;
  elec.R[1][2] = 0.0;
  elec.createWalkers(1);

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  elec.createSK();

  elec.addTable(ions);
  elec.update();


  elec.allocateGPU(10000);
  elec.copyWalkersToGPU();
  elec.updateLists_GPU();

  CoulombPBCAB_CUDA cab(ions, elec);

  // Background charge term
  double consts = cab.evalConsts(elec);
  CHECK(consts == Approx(0.0267892759 * 4));

  double val = cab.evaluate(elec);
  //cout << "CUDA BCC H val = " << val << std::endl;
  CHECK(val == Approx(-2.219665062 + 0.0267892759 * 4)); // not validated

  // actual code path use addEnergy
  std::vector<double> local_energy(1);
  cab.addEnergy(elec, local_energy);
  //cout << "addEnergy = " << local_energy[0] << std::endl;

  // not validated, which means it is just copied from the QMCPACK output, and has not been
  // computed independently.  It should be the same as results of the call to 'evaluate' above,
  //  which use the CPU code path.
  CHECK(local_energy[0] == Approx(-2.219665062 + 0.0267892759 * 4));
}

TEST_CASE("Coulomb PBC A-A CUDA BCC H", "[hamiltonian][CUDA]")
{
  LRCoulombSingleton::CoulombHandler = 0;
  LRCoulombSingleton::this_lr_type   = LRCoulombSingleton::ESLER;

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.R.diagonal(3.77945227);
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  MCWalkerConfiguration ions(simulation_cell);
  //ParticleSet elec;

  ions.setName("ion");
  std::vector<int> agroup(1);
  agroup[0] = 2;
  ions.create(agroup);
  //ions.create({2});
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.R[1][0] = 1.88972614;
  ions.R[1][1] = 1.88972614;
  ions.R[1][2] = 1.88972614;
  ions.createWalkers(1);

  SpeciesSet& ion_species       = ions.getSpeciesSet();
  int pIdx                      = ion_species.addSpecies("H");
  int pChargeIdx                = ion_species.addAttribute("charge");
  ion_species(pChargeIdx, pIdx) = 1;
  ions.createSK();
  ions.update();


  ions.allocateGPU(10000);
  ions.copyWalkersToGPU();
  ions.updateLists_GPU();


  CoulombPBCAA_CUDA caa(ions, false);

  // Background charge term
  double consts = caa.evalConsts();
  CHECK(consts == Approx(-1.6752294515)); // not validated

  double val = caa.evaluate(ions);
  //cout << "CUDA BCC A-A ion H val = " << val << std::endl;
  CHECK(val == Approx(-0.9628996199)); // not validated

  // actual code path uses addEnergy
  std::vector<double> local_energy(1);
  caa.addEnergy(ions, local_energy);
  //cout << "addEnergy = " << local_energy[0] << std::endl;
  CHECK(local_energy[0] == Approx(-0.9628996199)); // not validated
}


} // namespace qmcplusplus
