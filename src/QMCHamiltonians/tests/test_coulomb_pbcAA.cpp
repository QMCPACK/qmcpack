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
#include "QMCHamiltonians/CoulombPBCAA.h"
#include <ResourceCollection.h>
#include "QMCWaveFunctions/TrialWaveFunction.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
TEST_CASE("Coulomb PBC A-A", "[hamiltonian]")
{
  const double vmad_sc               = -1.4186487397403098;
  LRCoulombSingleton::CoulombHandler = 0;

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.R.diagonal(1.0);
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet ions(simulation_cell);

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


  CoulombPBCAA caa(ions, false, false, false);

  // Background charge term
  double consts = caa.evalConsts();
  CHECK(consts == Approx(-3.1151210154));

  double val = caa.evaluate(ions);
  //std::cout << "val = " << val << std::endl;
  CHECK(val == Approx(vmad_sc));

  // supercell Madelung energy
  val = caa.MC0;
  CHECK(val == Approx(vmad_sc));
}

TEST_CASE("Coulomb PBC A-A BCC H", "[hamiltonian]")
{
  const double alat                  = 3.77945227;
  const double vmad_sc               = -1.4186487397403098 / alat;
  LRCoulombSingleton::CoulombHandler = 0;

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.R.diagonal(alat);
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet ions(simulation_cell);
  ParticleSet elec(simulation_cell);

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


  CoulombPBCAA caa(ions, false, false, false);

  // Background charge term
  double consts = caa.evalConsts();
  CHECK(consts == Approx(-1.675229452)); // not validated

  double val = caa.evaluate(elec);
  CHECK(val == Approx(-0.9628996199)); // not validated

  // supercell Madelung energy
  val = caa.MC0;
  CHECK(val == Approx(vmad_sc));
}

TEST_CASE("Coulomb PBC A-A elec", "[hamiltonian]")
{
  LRCoulombSingleton::CoulombHandler = 0;

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.R.diagonal(1.0);
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet elec(simulation_cell);

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
  elec.update();


  CoulombPBCAA caa(elec, true, false, false);

  // Self-energy correction, no background charge for e-e interaction
  double consts = caa.evalConsts();
  CHECK(consts == Approx(-3.1151210154));

  double val = caa.evaluate(elec);
  CHECK(val == Approx(-1.418648723)); // not validated
}

TEST_CASE("Coulomb PBC A-A BCC", "[hamiltonian]")
{
  const double alat                  = 1.0;
  const double vmad_bcc              = -1.819616724754322 / alat;
  LRCoulombSingleton::CoulombHandler = 0;

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.R         = 0.5 * alat;
  lattice.R(0, 0)   = -0.5 * alat;
  lattice.R(1, 1)   = -0.5 * alat;
  lattice.R(2, 2)   = -0.5 * alat;
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet elec(simulation_cell);

  elec.setName("elec");
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


  CoulombPBCAA caa(elec, true, false, false);

  double val = caa.evaluate(elec);
  CHECK(val == Approx(vmad_bcc));

  val = caa.MC0;
  CHECK(val == Approx(vmad_bcc));
}

void test_CoulombPBCAA_3p(DynamicCoordinateKind kind)
{
  const double alat                  = 1.0;
  const double vmad_bcc              = -1.819616724754322 / alat;
  LRCoulombSingleton::CoulombHandler = 0;

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.R         = 0.5 * alat;
  lattice.R(0, 0)   = -0.5 * alat;
  lattice.R(1, 1)   = -0.5 * alat;
  lattice.R(2, 2)   = -0.5 * alat;
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet elec(simulation_cell, kind);

  elec.setName("elec");
  elec.create({1, 2});
  elec.R[0] = {0.0, 0.0, 0.0};
  elec.R[1] = {0.1, 0.2, 0.3};
  elec.R[2] = {0.3, 0.1, 0.2};

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int dnIdx                  = tspecies.addSpecies("d");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(chargeIdx, dnIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;
  tspecies(massIdx, dnIdx)   = 1.0;

  elec.createSK();

  CoulombPBCAA caa(elec, true, false, kind == DynamicCoordinateKind::DC_POS_OFFLOAD);

  ParticleSet elec_clone(elec);
  CoulombPBCAA caa_clone(caa);

  elec_clone.R[2] = {0.2, 0.3, 0.0};

  // testing batched interfaces
  ResourceCollection pset_res("test_pset_res");
  elec.createResource(pset_res);

  // testing batched interfaces
  RefVectorWithLeader<ParticleSet> p_ref_list(elec, {elec, elec_clone});
  RefVectorWithLeader<OperatorBase> caa_ref_list(caa, {caa, caa_clone});

  // dummy psi
  TrialWaveFunction psi, psi_clone;
  RefVectorWithLeader<TrialWaveFunction> psi_ref_list(psi, {psi, psi_clone});

  ResourceCollectionTeamLock<ParticleSet> mw_pset_lock(pset_res, p_ref_list);

  ParticleSet::mw_update(p_ref_list);
  caa.mw_evaluate(caa_ref_list, psi_ref_list, p_ref_list);

  CHECK(caa.getValue() == Approx(-5.4954533536));
  CHECK(caa_clone.getValue() == Approx(-6.329373489));

  CHECK(caa.MC0 == Approx(vmad_bcc));
  CHECK(caa_clone.MC0 == Approx(vmad_bcc));
}

TEST_CASE("Coulomb PBC A-A BCC 3 particles", "[hamiltonian]")
{
  test_CoulombPBCAA_3p(DynamicCoordinateKind::DC_POS);
  test_CoulombPBCAA_3p(DynamicCoordinateKind::DC_POS_OFFLOAD);
}
} // namespace qmcplusplus
