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
#include "Lattice/ParticleBConds.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "Particle/SymmetricDistanceTableData.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCApp/ParticleSetPool.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/ConstantOrbital.h"
#include "QMCWaveFunctions/LinearOrbital.h"
#include "QMCHamiltonians/BareKineticEnergy.h"
#include "Estimators/EstimatorManagerBase.h"
#include "Estimators/TraceManager.h"
#include "QMCDrivers/DMC/DMCUpdatePbyP.h"


#include <stdio.h>
#include <string>


using std::string;

namespace qmcplusplus
{

TEST_CASE("DMC Particle-by-Particle advanceWalkers ConstantOrbital", "[drivers][dmc]")
{

  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  ParticleSet ions;
  MCWalkerConfiguration elec;

  ions.setName("ion");
  ions.create(1);
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;

  elec.setName("elec");
  elec.setBoundBox(false);
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

  SpeciesSet &tspecies =  elec.getSpeciesSet();
  int upIdx = tspecies.addSpecies("u");
  int downIdx = tspecies.addSpecies("d");
  int chargeIdx = tspecies.addAttribute("charge");
  int massIdx = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(chargeIdx, downIdx) = -1;
  tspecies(massIdx, upIdx) = 1.0;
  tspecies(massIdx, downIdx) = 1.0;

  elec.addTable(ions,DT_AOS);
  elec.update();


  TrialWaveFunction psi(c);
  ConstantOrbital *orb = new ConstantOrbital;
  psi.addOrbital(orb, "Constant");
  psi.registerData(elec, elec.WalkerList[0]->DataSet);
  elec.WalkerList[0]->DataSet.allocate();

  FakeRandom rg;

  QMCHamiltonian h;
  h.addOperator(new BareKineticEnergy<double>(elec),"Kinetic");
  h.addObservables(elec); // get double free error on 'h.Observables' w/o this

  elec.resetWalkerProperty(); // get memory corruption w/o this

  DMCUpdatePbyPWithRejectionFast dmc(elec, psi, h, rg);
  EstimatorManagerBase EM;
  double tau = 0.1;
  SimpleFixedNodeBranch branch(tau, 1);
  TraceManager TM;
  dmc.resetRun(&branch, &EM, &TM);
  dmc.startBlock(1);

  DMCUpdatePbyPWithRejectionFast::WalkerIter_t begin = elec.begin();
  DMCUpdatePbyPWithRejectionFast::WalkerIter_t end = elec.end();
  dmc.advanceWalkers(begin, end, true);

  // With the constant wavefunction, no moves should be rejected
  REQUIRE(dmc.nReject == 0);
  REQUIRE(dmc.nAccept == 2);

  // Each electron moved sqrt(tau)*gaussian_rng()
  //  See ParticleBase/tests/test_random_seq.cpp for the gaussian random numbers
  REQUIRE(elec.R[0][0] == Approx(0.6276702589209545));
  REQUIRE(elec.R[0][1] == Approx(0.0));
  REQUIRE(elec.R[0][2] == Approx(-0.3723297410790455));

  REQUIRE(elec.R[1][0] == Approx(0.0));
  REQUIRE(elec.R[1][1] == Approx(-0.3723297410790455));
  REQUIRE(elec.R[1][2] == Approx(1.0));


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

  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  ParticleSet ions;
  MCWalkerConfiguration elec;

  ions.setName("ion");
  ions.create(1);
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;

  elec.setName("elec");
  elec.setBoundBox(false);
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

  SpeciesSet &tspecies =  elec.getSpeciesSet();
  int upIdx = tspecies.addSpecies("u");
  int downIdx = tspecies.addSpecies("d");
  int chargeIdx = tspecies.addAttribute("charge");
  int massIdx = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(chargeIdx, downIdx) = -1;
  tspecies(massIdx, upIdx) = 1.0;
  tspecies(massIdx, downIdx) = 1.0;

  elec.addTable(ions,DT_AOS);
  elec.update();


  TrialWaveFunction psi(c);
  LinearOrbital *orb = new LinearOrbital;
  psi.addOrbital(orb, "Linear");
  psi.registerData(elec, elec.WalkerList[0]->DataSet);
  elec.WalkerList[0]->DataSet.allocate();

  FakeRandom rg;

  QMCHamiltonian h;
  h.addOperator(new BareKineticEnergy<double>(elec),"Kinetic");
  h.addObservables(elec); // get double free error on 'h.Observables' w/o this

  elec.resetWalkerProperty(); // get memory corruption w/o this

  DMCUpdatePbyPWithRejectionFast dmc(elec, psi, h, rg);
  EstimatorManagerBase EM;
  double tau = 0.1;
  SimpleFixedNodeBranch branch(tau, 1);
  TraceManager TM;
  dmc.resetRun(&branch, &EM, &TM);
  dmc.startBlock(1);

  DMCUpdatePbyPWithRejectionFast::WalkerIter_t begin = elec.begin();
  DMCUpdatePbyPWithRejectionFast::WalkerIter_t end = elec.end();
  dmc.advanceWalkers(begin, end, true);

  REQUIRE(dmc.nReject == 0);
  REQUIRE(dmc.nAccept == 2);

  // Each electron moved sqrt(tau)*gaussian_rng()
  //  See ParticleBase/tests/test_random_seq.cpp for the gaussian random numbers
  //  See DMC_propagator notebook for computation of these values
  REQUIRE(elec.R[0][0] == Approx(0.695481606677082));
  REQUIRE(elec.R[0][1] == Approx(0.135622695565971));
  REQUIRE(elec.R[0][2] == Approx(-0.168895697756948));

  REQUIRE(elec.R[1][0] == Approx(0.0678113477829853));
  REQUIRE(elec.R[1][1] == Approx(-0.236707045539933));
  REQUIRE(elec.R[1][2] == Approx(1.20343404334896));

}
}

