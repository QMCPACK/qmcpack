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
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/ConstantOrbital.h"
#include "QMCHamiltonians/BareKineticEnergy.h"
#include "Estimators/EstimatorManagerBase.h"
#include "Estimators/TraceManager.h"
#include "QMCDrivers/VMC/VMCUpdatePbyP.h"
#include "QMCDrivers/GreenFunctionModifiers/DriftModifierUNR.h"


#include <stdio.h>
#include <string>


using std::string;

namespace qmcplusplus
{
TEST_CASE("VMC Particle-by-Particle advanceWalkers", "[drivers][vmc]")
{
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell);
  MCWalkerConfiguration elec(simulation_cell);

  ions.setName("ion");
  ions.create(1);
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;

  elec.setName("elec");
  std::vector<int> agroup(1);
  agroup[0] = 2;
  elec.create(agroup);
  elec.R[0][0] = 1.0;
  elec.R[0][1] = 0.0;
  elec.R[0][2] = 0.0;
  elec.R[1][0] = 0.0;
  elec.R[1][1] = 0.0;
  elec.R[1][2] = 1.0;
  elec.createWalkers(1);

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  elec.addTable(ions);
  elec.update();


  TrialWaveFunction psi;
  psi.addComponent(std::make_unique<ConstantOrbital>());
  psi.registerData(elec, elec.WalkerList[0]->DataSet);
  elec.WalkerList[0]->DataSet.allocate();

  FakeRandom rg;

  QMCHamiltonian h;
  h.addOperator(std::make_unique<BareKineticEnergy>(elec), "Kinetic");
  h.addObservables(elec); // get double free error on 'h.Observables' w/o this

  elec.resetWalkerProperty(); // get memory corruption w/o this

  VMCUpdatePbyP vmc(elec, psi, h, rg);
  EstimatorManagerBase EM;
  SimpleFixedNodeBranch branch(0.1, 1);
  TraceManager TM;
  DriftModifierUNR DM;
  vmc.resetRun(&branch, &EM, &TM, &DM);
  vmc.startBlock(1);

  VMCUpdatePbyP::WalkerIter_t begin = elec.begin();
  VMCUpdatePbyP::WalkerIter_t end   = elec.end();
  vmc.advanceWalkers(begin, end, true);

  // With the constant wavefunction, no moves should be rejected
  REQUIRE(vmc.nReject == 0);
  REQUIRE(vmc.nAccept == 2);

  // Each electron moved sqrt(tau)*gaussian_rng()
  //  See ParticleBase/tests/test_random_seq.cpp for the gaussian random numbers
  //  Values from diffuse.py
  REQUIRE(elec.R[0][0] == Approx(0.627670258894097));
  REQUIRE(elec.R[0][1] == Approx(0.0));
  REQUIRE(elec.R[0][2] == Approx(-0.372329741105903));

  REQUIRE(elec.R[1][0] == Approx(0.0));
  REQUIRE(elec.R[1][1] == Approx(-0.372329741105903));
  REQUIRE(elec.R[1][2] == Approx(1.0));
}
} // namespace qmcplusplus
