//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "EDenEstimatorManagerIntegrationTest.h"
#include "EnergyDensityTest.h"
#include "EstimatorManagerInput.h"
#include "EstimatorManagerInputTest.h"
#include "EstimatorInputDelegates.h"
#include "EstimatorManagerNew.h"
#include "EstimatorManagerCrowd.h"

namespace qmcplusplus
{
namespace testing
{

constexpr bool generate_test_data = false;
using MCPWalker                   = EDenEstimatorManagerIntegrationTest::MCPWalker;

EDenEstimatorManagerIntegrationTest::EDenEstimatorManagerIntegrationTest(Communicate* comm, int num_walkers)
{
  eden_test_ = std::make_unique<EnergyDensityTest>(comm, num_walkers, &testing::makeGoldWalkerElementsWithEEEIPS,
                                                   generate_test_data);
  auto doc   = testing::createEstimatorManagerEnergyDenistyInputXML();
  EstimatorManagerInput emi(doc.getRoot());
  auto& gold_elem = eden_test_->getGoldElements();
  auto ham_list   = eden_test_->getHamList();

  auto ham_lock = ResourceCollectionTeamLock<QMCHamiltonian>(eden_test_->getHamRes(), ham_list);

  auto pset_list = eden_test_->getPSetList();
  auto pset_lock = ResourceCollectionTeamLock<ParticleSet>(eden_test_->getPSetRes(), pset_list);

  auto twf_list = eden_test_->getTwfList();
  auto twf_lock = ResourceCollectionTeamLock<TrialWaveFunction>(eden_test_->getTwfRes(), twf_list);

  emn_ = std::make_unique<EstimatorManagerNew>(gold_elem.ham, comm);

  emn_->startDriverRun();
  emn_->constructEstimators(std::move(emi), gold_elem.pset_elec, gold_elem.twf, gold_elem.ham,
                            gold_elem.particle_pool.getPool());
  emc_ = std::make_unique<EstimatorManagerCrowd>(*emn_);
  emc_->registerListeners(ham_list);

  ParticleSet::mw_update(pset_list);
  TrialWaveFunction::mw_evaluateLog(twf_list, pset_list);
  QMCHamiltonian::mw_evaluate(ham_list, twf_list, pset_list);

  rng_.init(101);
}

void EDenEstimatorManagerIntegrationTest::updateAndEvaluate() {
  auto ham_list   = eden_test_->getHamList();

  auto ham_lock = ResourceCollectionTeamLock<QMCHamiltonian>(eden_test_->getHamRes(), ham_list);

  auto pset_list = eden_test_->getPSetList();
  auto pset_lock = ResourceCollectionTeamLock<ParticleSet>(eden_test_->getPSetRes(), pset_list);

  auto twf_list = eden_test_->getTwfList();
  auto twf_lock = ResourceCollectionTeamLock<TrialWaveFunction>(eden_test_->getTwfRes(), twf_list);

  ParticleSet::mw_update(pset_list);
  TrialWaveFunction::mw_evaluateLog(twf_list, pset_list);
  QMCHamiltonian::mw_evaluate(ham_list, twf_list, pset_list);
}
  
RefVectorWithLeader<MCPWalker> EDenEstimatorManagerIntegrationTest::getWalkerList()
{
  return eden_test_->getWalkerList();
}

RefVectorWithLeader<ParticleSet> EDenEstimatorManagerIntegrationTest::getPSetList()
{
  return eden_test_->getPSetList();
}

RefVectorWithLeader<QMCHamiltonian> EDenEstimatorManagerIntegrationTest::getHamList()
{
  return eden_test_->getHamList();
}

RefVectorWithLeader<TrialWaveFunction> EDenEstimatorManagerIntegrationTest::getTwfList()
{
  return eden_test_->getTwfList();
}

MockGoldWalkerElements& EDenEstimatorManagerIntegrationTest::getGoldElements() { return eden_test_->getGoldElements(); }

EstimatorManagerNew& EDenEstimatorManagerIntegrationTest::getEmn() { return *emn_; }

EstimatorManagerCrowd& EDenEstimatorManagerIntegrationTest::getEmc() { return *emc_; }

StdRandom<double>& EDenEstimatorManagerIntegrationTest::getRng() { return rng_; }
} // namespace testing
} // namespace qmcplusplus
