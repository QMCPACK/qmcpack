//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//
// File created by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//////////////////////////////////////////////////////////////////////////////////////
#include "catch.hpp"

#include "Particle/ParticleSet.h"
#include "QMCHamiltonians/CoulombPBCAA.h"
#include "LongRange/EwaldHandler2D.h"

namespace qmcplusplus
{

TEST_CASE("Coulomb PBC A-A Ewald2D square", "[hamiltonian]")
{
  LRCoulombSingleton::CoulombHandler = 0;
  LRCoulombSingleton::this_lr_type = LRCoulombSingleton::STRICT2D;
  const double vmad_sq = -1.95013246;
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.ndim = 2;
  lattice.R.diagonal(1.0);
  lattice.LR_dim_cutoff = 30.0;
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet elec(simulation_cell);
  elec.setName("e");
  elec.create({1});
  elec.R[0] = {0.0, 0.0, 0.0};

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  elec.createSK();
  elec.update();

  CoulombPBCAA caa = CoulombPBCAA(elec, true, false, false);

  double val = caa.evaluate(elec);
  CHECK(val == Approx(vmad_sq));
}

TEST_CASE("Coulomb PBC A-A Ewald2D body center", "[hamiltonian]")
{
  LRCoulombSingleton::CoulombHandler = 0; // !!!! crucial if not first test
  LRCoulombSingleton::this_lr_type = LRCoulombSingleton::STRICT2D;
  const double vmad_bc = -2.7579038;

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.ndim = 2;
  lattice.R = 0.0;
  lattice.R.diagonal(1.0);
  lattice.LR_dim_cutoff = 30.0;
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet elec(simulation_cell);
  elec.setName("e");
  const int npart = 2;
  elec.create({npart});
  // initialize fractional coordinates
  elec.R[0] = {0.0, 0.0, 0.0};
  elec.R[1] = {0.5, 0.5, 0.0};
  // convert to Cartesian coordinates
  for (int i=0;i<npart;i++)
    elec.R[i] = dot(elec.R[i], lattice.R);

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  elec.createSK();
  elec.addTable(elec); // !!!! crucial with more than 1 particle
  elec.update();

  CoulombPBCAA caa = CoulombPBCAA(elec, true, false, false);

  double val = caa.evaluate(elec);
  CHECK(val/npart == Approx(vmad_bc));
}

TEST_CASE("Coulomb PBC A-A Ewald2D triangular", "[hamiltonian]")
{
  LRCoulombSingleton::CoulombHandler = 0; // !!!! crucial if not first test
  LRCoulombSingleton::this_lr_type = LRCoulombSingleton::STRICT2D;
  const double vmad_tri = -1.1061025865191676;
  const double alat = std::sqrt(2.0*M_PI/std::sqrt(3));

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.ndim = 2;
  lattice.R = 0.0;
  lattice.R(0, 0) = alat;
  lattice.R(1, 0) = -1.0/2*alat;
  lattice.R(1, 1) = std::sqrt(3)/2*alat;
  lattice.R(2, 2) = 2*alat;
  lattice.LR_dim_cutoff = 30.0;
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet elec(simulation_cell);
  elec.setName("e");
  elec.create({1});
  elec.R[0] = {0.0, 0.0, 0.0};

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  elec.createSK();
  elec.update();

  CoulombPBCAA caa = CoulombPBCAA(elec, true, false, false);

  double val = caa.evaluate(elec);
  CHECK(val == Approx(vmad_tri));
}

TEST_CASE("Coulomb PBC A-A Ewald2D honeycomb", "[hamiltonian]")
{
  LRCoulombSingleton::CoulombHandler = 0; // !!!! crucial if not first test
  LRCoulombSingleton::this_lr_type = LRCoulombSingleton::STRICT2D;
  const double vmad_hon = -1.510964233;
  const double alat = std::sqrt(2.0*M_PI/std::sqrt(3));

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.ndim = 2;
  lattice.R = 0.0;
  lattice.R(0, 0) = alat;
  lattice.R(1, 0) = -1.0/2*alat;
  lattice.R(1, 1) = std::sqrt(3)/2*alat;
  lattice.R(2, 2) = 2*alat;
  lattice.LR_dim_cutoff = 30.0;
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet elec(simulation_cell);
  elec.setName("e");
  const int npart = 2;
  elec.create({npart});
  // initialize fractional coordinates
  elec.R[0] = {0.0, 0.0, 0.0};
  elec.R[1] = {2./3, 1./3, 0.0};
  // convert to Cartesian coordinates
  for (int i=0;i<npart;i++)
    elec.R[i] = dot(elec.R[i], lattice.R);

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  elec.createSK();
  elec.addTable(elec); // !!!! crucial with more than 1 particle
  elec.update();

  CoulombPBCAA caa = CoulombPBCAA(elec, true, false, false);

  double val = caa.evaluate(elec);
  CHECK(val/npart == Approx(vmad_hon));
}

TEST_CASE("Coulomb PBC A-A Ewald2D tri. in rect.", "[hamiltonian]")
{
  LRCoulombSingleton::CoulombHandler = 0; // !!!! crucial if not first test
  LRCoulombSingleton::this_lr_type = LRCoulombSingleton::STRICT2D;
  const double vmad_tri = -1.1061025865191676;
  const double alat = std::sqrt(2.0*M_PI/std::sqrt(3));

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.ndim = 2;
  lattice.R = 0.0;
  lattice.R(0, 0) = alat;
  lattice.R(1, 1) = std::sqrt(3)*alat;
  lattice.R(2, 2) = 4.0;
  lattice.LR_dim_cutoff = 30.0;
  lattice.reset();
  lattice.print(app_log());

  const SimulationCell simulation_cell(lattice);
  ParticleSet elec(simulation_cell);
  elec.setName("e");
  elec.create({2});
  elec.R[0] = {0.0, 0.0, 0.0};
  elec.R[1] = {alat/2, std::sqrt(3.0)*alat/2, 0.0};

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  elec.createSK();
  elec.addTable(elec); // !!!! crucial with more than 1 particle
  elec.update();
  app_log() << elec.R << std::endl;

  CoulombPBCAA caa = CoulombPBCAA(elec, true, false, false);

  double val = caa.evaluate(elec)/elec.getTotalNum();
  CHECK(val == Approx(vmad_tri));
}
} // qmcplusplus
