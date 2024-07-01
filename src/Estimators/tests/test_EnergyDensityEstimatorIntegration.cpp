//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "EstimatorManagerNew.h"
#include "EnergyDensityTest.h"
#include "EnergyDensityEstimator.h"
#include "EstimatorManagerInputTest.h"
#include "EstimatorManagerInput.h"
#include "EstimatorInputDelegates.h"
#include "EstimatorManagerCrowd.h"
#include "Utilities/ProjectData.h"
#include "GenerateRandomParticleSets.h"

constexpr bool generate_test_data = true;

namespace qmcplusplus
{

TEST_CASE("EnergyDensityEstimatorIntegration", "[estimators]")
{
  Communicate* comm = OHMMS::Controller;

  testing::EnergyDensityTest eden_test(comm, 4 /*num_walkers*/, generate_test_data);

  auto doc = testing::createEstimatorManagerEnergyDenistyInputXML();
  EstimatorManagerInput emi(doc.getRoot());
  auto& gold_elem = eden_test.getGoldElements();
  EstimatorManagerNew emn(comm, std::move(emi), gold_elem.ham, gold_elem.pset_elec, gold_elem.particle_pool.getPool(), gold_elem.twf);

  auto pset_list = eden_test.getPSetList();
  auto pset_lock = eden_test.makeTeamLock(eden_test.getPSetRes(), pset_list);
  pset_list[0].L[0] = 1.0;
  pset_list[1].L[1] = 1.0;
  pset_list[2].L[2] = 1.0;
  pset_list[3].L[3] = 1.0;

  EstimatorManagerCrowd emc(emn);  

  auto ham_list = eden_test.getHamList();
  auto ham_lock = eden_test.makeTeamLock(eden_test.getHamRes(), ham_list);
  emc.registerListeners(ham_list);

  auto twf_list = eden_test.getTwfList();
  auto twf_lock = eden_test.makeTeamLock(eden_test.getTwfRes(), twf_list);

  ParticleSet::mw_update(pset_list);
  TrialWaveFunction::mw_evaluateLog(twf_list, pset_list);
  QMCHamiltonian::mw_evaluate(ham_list, twf_list, pset_list);

  StdRandom<double> rng;
  rng.init(101);

  auto walker_list = eden_test.getWalkerList();
  emc.accumulate(walker_list, pset_list, twf_list, ham_list, rng);
}

}
