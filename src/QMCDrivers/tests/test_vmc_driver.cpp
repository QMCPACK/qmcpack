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
#include "Lattice/ParticleBConds.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "Particle/SymmetricDistanceTableData.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCApp/ParticleSetPool.h"
#include "QMCApp/HamiltonianPool.h"
#include "QMCApp/WaveFunctionPool.h"
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
  Communicate* c;
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

  SpeciesSet& tspecies         = elec.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int downIdx                  = tspecies.addSpecies("d");
  int chargeIdx                = tspecies.addAttribute("charge");
  int massIdx                  = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;
  tspecies(massIdx, upIdx)     = 1.0;
  tspecies(massIdx, downIdx)   = 1.0;

#ifdef ENABLE_SOA
  elec.addTable(ions, DT_SOA);
#else
  elec.addTable(ions, DT_AOS);
#endif
  elec.update();

  CloneManager::clear_for_unit_tests();

  TrialWaveFunction psi = TrialWaveFunction(c);
  ConstantOrbital* orb  = new ConstantOrbital;
  psi.addComponent(orb, "Constant");
  psi.registerData(elec, elec.WalkerList[0]->DataSet);
  elec.WalkerList[0]->DataSet.allocate();

  FakeRandom rg;

  QMCHamiltonian h;
  h.addOperator(new BareKineticEnergy<double>(elec), "Kinetic");
  h.addObservables(elec); // get double free error on 'h.Observables' w/o this

  elec.resetWalkerProperty(); // get memory corruption w/o this

  HamiltonianPool hpool(c);

  WaveFunctionPool wpool(c);

  //EstimatorManagerBase emb(c);


  VMC vmc_omp(elec, psi, h, wpool, c);

  const char* vmc_input = "<qmc method=\"vmc\" move=\"pbyp\"> \
   <parameter name=\"substeps\">1</parameter> \
   <parameter name=\"steps\">1</parameter> \
   <parameter name=\"blocks\">1</parameter> \
   <parameter name=\"timestep\">0.1</parameter> \
   <parameter name=\"usedrift\">no</parameter> \
  </qmc> \
  ";
  Libxml2Document* doc  = new Libxml2Document;
  bool okay             = doc->parseFromString(vmc_input);
  REQUIRE(okay);
  xmlNodePtr root = doc->getRoot();

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
} // namespace qmcplusplus
