//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"


#include "Utilities/RandomGenerator.h"
#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/ParticleSetPool.h"
#include "QMCHamiltonians/HamiltonianPool.h"
#include "QMCWaveFunctions/WaveFunctionPool.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/ConstantOrbital.h"
#include "QMCHamiltonians/BareKineticEnergy.h"
#include "Estimators/EstimatorManagerBase.h"
#include "Estimators/TraceManager.h"
#include "QMCDrivers/VMC/VMC.h"


#include <stdio.h>
#include <string>


using std::string;

namespace qmcplusplus
{
TEST_CASE("VMC", "[drivers][vmc]")
{
  Communicate* c = OHMMS::Controller;
  c->setName("test");
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell);
  MCWalkerConfiguration elec(simulation_cell);

  ions.setName("ion");
  ions.create(1);
  ions.R[0] = {0.0, 0.0, 0.0};

  elec.setName("elec");
  std::vector<int> agroup(1, 2);
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

  CloneManager::clearClones();

  TrialWaveFunction psi;
  psi.addComponent(std::make_unique<ConstantOrbital>());
  psi.registerData(elec, elec.WalkerList[0]->DataSet);
  elec.WalkerList[0]->DataSet.allocate();

  FakeRandom rg;

  QMCHamiltonian h;
  h.addOperator(std::make_unique<BareKineticEnergy>(elec), "Kinetic");
  h.addObservables(elec); // get double free error on 'h.Observables' w/o this

  elec.resetWalkerProperty(); // get memory corruption w/o this

  VMC vmc_omp(elec, psi, h, c, false);

  const char* vmc_input = "<qmc method=\"vmc\" move=\"pbyp\"> \
   <parameter name=\"substeps\">1</parameter> \
   <parameter name=\"steps\">1</parameter> \
   <parameter name=\"blocks\">1</parameter> \
   <parameter name=\"timestep\">0.1</parameter> \
   <parameter name=\"usedrift\">no</parameter> \
  </qmc> \
  ";
  Libxml2Document doc;
  bool okay = doc.parseFromString(vmc_input);
  REQUIRE(okay);
  xmlNodePtr root = doc.getRoot();

  vmc_omp.process(root); // need to call 'process' for QMCDriver, which in turn calls 'put'

  vmc_omp.run();

  // With the constant wavefunction, no moves should be rejected
  double ar = vmc_omp.acceptRatio();
  REQUIRE(ar == Approx(1.0));

  // Each electron moved sqrt(tau)*gaussian_rng()
  //  See Particle>Base/tests/test_random_seq.cpp for the gaussian random numbers
  //  Values from diffuse.py for moving one step

  REQUIRE(elec[0]->R[0][0] == Approx(0.627670258894097));
  REQUIRE(elec.R[0][1] == Approx(0.0));
  REQUIRE(elec.R[0][2] == Approx(-0.372329741105903));

  REQUIRE(elec.R[1][0] == Approx(0.0));
  REQUIRE(elec.R[1][1] == Approx(-0.372329741105903));
  REQUIRE(elec.R[1][2] == Approx(1.0));
}

TEST_CASE("SOVMC", "[drivers][vmc]")
{
  Communicate* c = OHMMS::Controller;
  c->setName("test");
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell);
  MCWalkerConfiguration elec(simulation_cell);

  ions.setName("ion");
  ions.create(1);
  ions.R[0] = {0.0, 0.0, 0.0};

  elec.setName("elec");
  std::vector<int> agroup(1, 1);
  elec.create(agroup);
  elec.R[0]     = {1.0, 0.0, 0.0};
  elec.spins[0] = 0.0;
  elec.setSpinor(true);
  elec.createWalkers(1);

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  elec.addTable(ions);
  elec.update();

  CloneManager::clearClones();

  TrialWaveFunction psi;
  psi.addComponent(std::make_unique<ConstantOrbital>());
  psi.registerData(elec, elec.WalkerList[0]->DataSet);
  elec.WalkerList[0]->DataSet.allocate();

  FakeRandom rg;

  QMCHamiltonian h;
  h.addOperator(std::make_unique<BareKineticEnergy>(elec), "Kinetic");
  h.addObservables(elec); // get double free error on 'h.Observables' w/o this

  elec.resetWalkerProperty(); // get memory corruption w/o this

  VMC vmc_omp(elec, psi, h, c, false);

  const char* vmc_input = "<qmc method=\"vmc\" move=\"pbyp\"> \
   <parameter name=\"substeps\">1</parameter> \
   <parameter name=\"steps\">1</parameter> \
   <parameter name=\"blocks\">1</parameter> \
   <parameter name=\"timestep\">0.1</parameter> \
   <parameter name=\"usedrift\">no</parameter> \
   <parameter name=\"SpinMass\">0.25</parameter> \
  </qmc> \
  ";
  Libxml2Document doc;
  bool okay = doc.parseFromString(vmc_input);
  REQUIRE(okay);
  xmlNodePtr root = doc.getRoot();

  vmc_omp.process(root); // need to call 'process' for QMCDriver, which in turn calls 'put'

  vmc_omp.run();

  // With the constant wavefunction, no moves should be rejected
  double ar = vmc_omp.acceptRatio();
  REQUIRE(ar == Approx(1.0));

  // Each electron moved sqrt(tau)*gaussian_rng()
  //  See Particle>Base/tests/test_random_seq.cpp for the gaussian random numbers
  //  Values from diffuse.py for moving one step

  REQUIRE(elec.R[0][0] == Approx(0.627670258894097));
  REQUIRE(elec.R[0][1] == Approx(0.0));
  REQUIRE(elec.R[0][2] == Approx(-0.372329741105903));

  REQUIRE(elec.spins[0] == Approx(-0.74465948215809097));

  //Now we're going to test that the step updated the walker variables.
  REQUIRE(elec.WalkerList[0]->R[0][0] == Approx(elec.R[0][0]));
  REQUIRE(elec.WalkerList[0]->R[0][1] == Approx(elec.R[0][1]));
  REQUIRE(elec.WalkerList[0]->R[0][2] == Approx(elec.R[0][2]));
  REQUIRE(elec.WalkerList[0]->spins[0] == Approx(elec.spins[0]));
}

TEST_CASE("SOVMC-alle", "[drivers][vmc]")
{
  Communicate* c = OHMMS::Controller;
  c->setName("test");
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell);
  MCWalkerConfiguration elec(simulation_cell);

  ions.setName("ion");
  ions.create(1);
  ions.R[0] = {0.0, 0.0, 0.0};

  elec.setName("elec");
  std::vector<int> agroup(1, 1);
  elec.create(agroup);
  elec.R[0]     = {1.0, 0.0, 0.0};
  elec.spins[0] = 0.0;
  elec.setSpinor(true);
  elec.createWalkers(1);

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  elec.addTable(ions);
  elec.update();

  CloneManager::clearClones();

  TrialWaveFunction psi;
  psi.addComponent(std::make_unique<ConstantOrbital>());
  psi.registerData(elec, elec.WalkerList[0]->DataSet);
  elec.WalkerList[0]->DataSet.allocate();

  FakeRandom rg;

  QMCHamiltonian h;
  h.addOperator(std::make_unique<BareKineticEnergy>(elec), "Kinetic");
  h.addObservables(elec); // get double free error on 'h.Observables' w/o this

  elec.resetWalkerProperty(); // get memory corruption w/o this

  VMC vmc_omp(elec, psi, h, c, false);

  const char* vmc_input = "<qmc method=\"vmc\" move=\"alle\"> \
   <parameter name=\"substeps\">1</parameter> \
   <parameter name=\"steps\">1</parameter> \
   <parameter name=\"blocks\">1</parameter> \
   <parameter name=\"timestep\">0.1</parameter> \
   <parameter name=\"usedrift\">no</parameter> \
   <parameter name=\"SpinMass\">0.25</parameter> \
  </qmc> \
  ";
  Libxml2Document doc;
  bool okay = doc.parseFromString(vmc_input);
  REQUIRE(okay);
  xmlNodePtr root = doc.getRoot();

  vmc_omp.process(root); // need to call 'process' for QMCDriver, which in turn calls 'put'

  vmc_omp.run();

  // With the constant wavefunction, no moves should be rejected
  double ar = vmc_omp.acceptRatio();
  REQUIRE(ar == Approx(1.0));

  // Each electron moved sqrt(tau)*gaussian_rng()
  //  See Particle>Base/tests/test_random_seq.cpp for the gaussian random numbers
  //  Values from diffuse.py for moving one step

  REQUIRE(elec.R[0][0] == Approx(0.627670258894097));
  REQUIRE(elec.R[0][1] == Approx(0.0));
  REQUIRE(elec.R[0][2] == Approx(-0.372329741105903));

  REQUIRE(elec.spins[0] == Approx(-0.74465948215809097));

  //Now we're going to test that the step updated the walker variables.
  REQUIRE(elec.WalkerList[0]->R[0][0] == Approx(elec.R[0][0]));
  REQUIRE(elec.WalkerList[0]->R[0][1] == Approx(elec.R[0][1]));
  REQUIRE(elec.WalkerList[0]->R[0][2] == Approx(elec.R[0][2]));
  REQUIRE(elec.WalkerList[0]->spins[0] == Approx(elec.spins[0]));
}
} // namespace qmcplusplus
