//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include <ResourceCollection.h>
#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "QMCHamiltonians/CoulombPBCAB.h"
#include "QMCHamiltonians/CoulombPBCAA.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "TestListenerFunction.h"
#include "Utilities/RuntimeOptions.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{

using QMCT = QMCTraits;
using Real = QMCT::RealType;

TEST_CASE("Coulomb PBC A-B", "[hamiltonian]")
{
  LRCoulombSingleton::CoulombHandler = 0;

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


  CoulombPBCAB cab(ions, elec);

  // Self energy plus Background charge term
  CHECK(cab.myConst == Approx(2 * 0.0506238028)); // not validated

  double val_ei = cab.evaluate(elec);
  CHECK(val_ei == Approx(-0.005314032183 + 2 * 0.0506238028)); // not validated

  CoulombPBCAA caa_elec(elec, true, false, false);
  CoulombPBCAA caa_ion(ions, false, false, false);
  double val_ee = caa_elec.evaluate(elec);
  double val_ii = caa_ion.evaluate(ions);
  double sum    = val_ee + val_ii + val_ei;
  CHECK(sum == Approx(-2.741363553)); // Can be validated via Ewald summation elsewhere
                                      // -2.74136517454081
}

TEST_CASE("Coulomb PBC A-B BCC H", "[hamiltonian]")
{
  LRCoulombSingleton::CoulombHandler = 0;

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

  CoulombPBCAB cab(ions, elec);

  // Background charge term
  double consts = cab.evalConsts(elec);
  CHECK(consts == Approx(0.0267892759 * 4)); // not validated


  double val_ei = cab.evaluate(elec);
  CHECK(val_ei == Approx(-2.219665062 + 0.0267892759 * 4)); // not validated


  CoulombPBCAA caa_elec(elec, false, false, false);
  CoulombPBCAA caa_ion(ions, false, false, false);
  double val_ee = caa_elec.evaluate(elec);
  double val_ii = caa_ion.evaluate(ions);
  double sum    = val_ee + val_ii + val_ei;
  CHECK(sum == Approx(-3.143491064)); // Can be validated via Ewald summation elsewhere
                                      // -3.14349127313640
}

TEST_CASE("CoulombAB::Listener", "[hamiltonian]")
{
  using testing::getParticularListener;

  LRCoulombSingleton::CoulombHandler = 0;

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.R.diagonal(3.77945227);
  lattice.reset();

  const SimulationCell simulation_cell(lattice);

  ParticleSet ions(simulation_cell);

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
  ions.turnOnPerParticleSK();

  ParticleSet ions2(ions);
  ions2.update();

  ParticleSet elec(simulation_cell);

  elec.setName("elec");
  elec.create({2});
  elec.R[0]                  = {0.0, 0.5, 0.0};
  elec.R[1]                  = {0.0, 0.5, 0.0};
  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  elec.createSK();
  elec.update();
  elec.turnOnPerParticleSK();

  ParticleSet elec2(elec);

  elec2.R[0] = {0.0, 0.5, 0.1};
  elec2.R[1] = {0.6, 0.05, -0.1};
  elec2.update();

  CoulombPBCAB cab(ions, elec, false);
  CoulombPBCAB cab2(ions2, elec2, false);
  RefVector<OperatorBase> cabs{cab, cab2};
  RefVectorWithLeader<OperatorBase> o_list(cab, cabs);
  // Self-energy correction, no background charge for e-e interaction
  double consts = cab.evalConsts(elec);
  CHECK(consts == Approx(0.0267892759 * 4));
  RefVector<ParticleSet> ptcls{elec, elec2};
  RefVectorWithLeader<ParticleSet> p_list(elec, ptcls);

  RefVector<ParticleSet> ion_ptcls{ions, ions2};
  RefVectorWithLeader<ParticleSet> ion_p_list(ions, ion_ptcls);

  RuntimeOptions runtime_options;
  TrialWaveFunction psi(runtime_options);
  TrialWaveFunction psi_clone(runtime_options);
  RefVectorWithLeader<TrialWaveFunction> twf_list(psi, {psi, psi_clone});

  Matrix<Real> local_pots(2);
  Matrix<Real> local_pots2(2);

  ResourceCollection cab_res("test_cab_res");
  cab.createResource(cab_res);
  ResourceCollectionTeamLock<OperatorBase> cab_lock(cab_res, o_list);

  ResourceCollection pset_res("test_pset_res");
  elec.createResource(pset_res);
  ResourceCollectionTeamLock<ParticleSet> pset_lock(pset_res, p_list);

  std::vector<ListenerVector<Real>> listeners;
  listeners.emplace_back("localpotential", getParticularListener(local_pots));
  listeners.emplace_back("localpotential", getParticularListener(local_pots2));

  Matrix<Real> ion_pots(2);
  Matrix<Real> ion_pots2(2);

  std::vector<ListenerVector<Real>> ion_listeners;
  ion_listeners.emplace_back("localpotential", getParticularListener(ion_pots));
  ion_listeners.emplace_back("localpotential", getParticularListener(ion_pots2));

  ParticleSet::mw_update(p_list);
  cab.mw_evaluatePerParticle(o_list, twf_list, p_list, listeners, ion_listeners);
  CHECK(cab.getValue() == Approx(-2.219665062 + 0.0267892759 * 4));
  CHECK(cab2.getValue() == Approx(-1.7222343352));
  // Check that the sum of the particle energies == the total
  Real elec_ion_sum = std::accumulate(local_pots.begin(), local_pots.begin() + local_pots.cols(), 0.0);
  CHECK(elec_ion_sum == Approx(-1.0562537047));
  CHECK(elec_ion_sum + std::accumulate(ion_pots.begin(), ion_pots.begin() + ion_pots.cols(), 0.0) ==
        Approx(-2.219665062 + 0.0267892759 * 4));
  elec_ion_sum = std::accumulate(local_pots[1], local_pots[1] + local_pots.cols(), 0.0);
  CHECK(elec_ion_sum == Approx(-0.8611171676));
  CHECK(elec_ion_sum + std::accumulate(ion_pots[1], ion_pots[1] + ion_pots.cols(), 0.0) == Approx(-1.7222343352));
}


} // namespace qmcplusplus
