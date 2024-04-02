//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Utilities/RandomGenerator.h"
#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/ParticleSetPool.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/ConstantOrbital.h"
#include "QMCWaveFunctions/LinearOrbital.h"
#include "QMCHamiltonians/BareKineticEnergy.h"
#include "Estimators/EstimatorManagerBase.h"
#include "Estimators/TraceManager.h"
#include "QMCDrivers/DMC/DMCUpdatePbyP.h"
#include "QMCDrivers/GreenFunctionModifiers/DriftModifierUNR.h"
#include "Utilities/RuntimeOptions.h"


#include <stdio.h>
#include <string>


using std::string;

namespace qmcplusplus
{
TEST_CASE("DMC Particle-by-Particle advanceWalkers ConstantOrbital", "[drivers][dmc]")
{
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell);
  MCWalkerConfiguration elec(simulation_cell);

  ions.setName("ion");
  ions.create({1});
  ions.R[0] = {0.0, 0.0, 0.0};
  elec.setName("elec");
  std::vector<int> agroup(1);
  agroup[0] = 2;
  elec.create(agroup);
  elec.R[0] = {1.0, 0.0, 0.0};
  elec.R[1] = {0.0, 0.0, 1.0};
  elec.createWalkers(1);

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  elec.addTable(ions);
  elec.update();

  RuntimeOptions runtime_options;
  TrialWaveFunction psi(runtime_options);
  auto orb_uptr = std::make_unique<ConstantOrbital>();
  auto orb      = orb_uptr.get();
  psi.addComponent(std::move(orb_uptr));
  psi.registerData(elec, elec[0]->DataSet);
  elec[0]->DataSet.allocate();

  FakeRandom rg;

  QMCHamiltonian h;
  h.addOperator(std::make_unique<BareKineticEnergy>(elec, psi), "Kinetic");
  h.addObservables(elec); // get double free error on 'h.Observables' w/o this

  elec.resetWalkerProperty(); // get memory corruption w/o this

  DMCUpdatePbyPWithRejectionFast dmc(elec, psi, h, rg);
  EstimatorManagerBase EM;
  double tau = 0.1;
  SimpleFixedNodeBranch branch(tau, 1);
  TraceManager TM;
  DriftModifierUNR DM;
  dmc.resetRun(&branch, &EM, &TM, &DM);
  dmc.startBlock(1);

  DMCUpdatePbyPWithRejectionFast::WalkerIter_t begin = elec.begin();
  DMCUpdatePbyPWithRejectionFast::WalkerIter_t end   = elec.end();
  dmc.advanceWalkers(begin, end, true);

  // With the constant wavefunction, no moves should be rejected
  REQUIRE(dmc.nReject == 0);
  REQUIRE(dmc.nAccept == 2);

  // Each electron moved sqrt(tau)*gaussian_rng()
  //  See ParticleBase/tests/test_random_seq.cpp for the gaussian random numbers
  CHECK(elec.R[0][0] == Approx(0.6276702589209545));
  CHECK(elec.R[0][1] == Approx(0.0));
  CHECK(elec.R[0][2] == Approx(-0.3723297410790455));

  CHECK(elec.R[1][0] == Approx(0.0));
  CHECK(elec.R[1][1] == Approx(-0.3723297410790455));
  CHECK(elec.R[1][2] == Approx(1.0));


  // Check rejection in case of node-crossing

  orb->FakeGradRatio = -1.0;

  dmc.advanceWalkers(begin, end, true);
#ifdef QMC_COMPLEX
  REQUIRE(dmc.nReject == 0);
  REQUIRE(dmc.nAccept == 4);
#else
  // Should be rejected
  REQUIRE(dmc.nReject == 2);
  REQUIRE(dmc.nAccept == 2);
#endif
}

TEST_CASE("DMC Particle-by-Particle advanceWalkers LinearOrbital", "[drivers][dmc]")

{
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell);
  MCWalkerConfiguration elec(simulation_cell);

  ions.setName("ion");
  ions.create({1});
  ions.R[0] = {0.0, 0.0, 0.0};
  elec.setName("elec");
  std::vector<int> agroup(1);
  agroup[0] = 2;
  elec.create(agroup);
  elec.R[0] = {1.0, 0.0, 0.0};
  elec.R[1] = {0.0, 0.0, 1.0};
  elec.createWalkers(1);

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  elec.addTable(ions);
  elec.update();

  RuntimeOptions runtime_options;
  TrialWaveFunction psi(runtime_options);
  psi.addComponent(std::make_unique<LinearOrbital>());
  psi.registerData(elec, elec[0]->DataSet);
  elec[0]->DataSet.allocate();

  FakeRandom rg;

  QMCHamiltonian h;
  h.addOperator(std::make_unique<BareKineticEnergy>(elec, psi), "Kinetic");
  h.addObservables(elec); // get double free error on 'h.Observables' w/o this

  elec.resetWalkerProperty(); // get memory corruption w/o this

  DMCUpdatePbyPWithRejectionFast dmc(elec, psi, h, rg);
  EstimatorManagerBase EM;
  double tau = 0.1;
  SimpleFixedNodeBranch branch(tau, 1);
  TraceManager TM;
  DriftModifierUNR DM;
  dmc.resetRun(&branch, &EM, &TM, &DM);
  dmc.startBlock(1);

  DMCUpdatePbyPWithRejectionFast::WalkerIter_t begin = elec.begin();
  DMCUpdatePbyPWithRejectionFast::WalkerIter_t end   = elec.end();
  dmc.advanceWalkers(begin, end, true);

  REQUIRE(dmc.nReject == 0);
  REQUIRE(dmc.nAccept == 2);

  // Each electron moved sqrt(tau)*gaussian_rng()
  //  See ParticleBase/tests/test_random_seq.cpp for the gaussian random numbers
  //  See DMC_propagator notebook for computation of these values
  CHECK(elec.R[0][0] == Approx(0.695481606677082));
  CHECK(elec.R[0][1] == Approx(0.135622695565971));
  CHECK(elec.R[0][2] == Approx(-0.168895697756948));

  CHECK(elec.R[1][0] == Approx(0.0678113477829853));
  CHECK(elec.R[1][1] == Approx(-0.236707045539933));
  CHECK(elec.R[1][2] == Approx(1.20343404334896));
}
} // namespace qmcplusplus
