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
#include "PerParticleHamiltonianLogger.h"
#include "EnergyDensityTest.h"
#include "EnergyDensityEstimator.h"
#include "EstimatorManagerInputTest.h"
#include "EstimatorManagerInput.h"
#include "EstimatorInputDelegates.h"
#include "EstimatorManagerCrowd.h"
#include "Utilities/ProjectData.h"
#include "GenerateRandomParticleSets.h"
#include "EstimatorManagerNewTest.h"

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
  auto ham_list   = eden_test.getHamList();

  auto ham_lock = ResourceCollectionTeamLock<QMCHamiltonian>(eden_test.getHamRes(), ham_list);

  auto pset_list    = eden_test.getPSetList();
  auto pset_lock    = ResourceCollectionTeamLock<ParticleSet>(eden_test.getPSetRes(), pset_list);
  pset_list[0].L[0] = 1.0;
  pset_list[1].L[1] = 1.0;
  pset_list[2].L[2] = 1.0;
  pset_list[3].L[3] = 1.0;

  auto twf_list = eden_test.getTwfList();
  auto twf_lock = ResourceCollectionTeamLock<TrialWaveFunction>(eden_test.getTwfRes(), twf_list);

  EstimatorManagerNew emn(gold_elem.ham, comm);
  emn.constructEstimators(std::move(emi), gold_elem.pset_elec, gold_elem.twf, gold_elem.ham,
                          gold_elem.particle_pool.getPool());
  EstimatorManagerCrowd emc(emn);
  emc.registerListeners(ham_list);

  ParticleSet::mw_update(pset_list);
  TrialWaveFunction::mw_evaluateLog(twf_list, pset_list);
  QMCHamiltonian::mw_evaluate(ham_list, twf_list, pset_list);

  StdRandom<double> rng;
  rng.init(101);

  auto walker_list = eden_test.getWalkerList();
  emc.accumulate(walker_list, pset_list, twf_list, ham_list, rng);

  std::vector<RefVector<OperatorEstBase>> crowd_operator_ests = {emc.get_operator_estimators()};
  emn.collectOperatorEstimators(crowd_operator_ests);

  testing::EstimatorManagerNewTestAccess emnta(emn);
  emnta.reduceOperatorEstimators();
  auto operator_ests = emnta.getOperatorEstimators();

  NEEnergyDensityEstimator& e_den_est = dynamic_cast<NEEnergyDensityEstimator&>(operator_ests[0].get());

  auto spacegrids = e_den_est.getSpaceGrids();

  decltype(spacegrids)::value_type::type& grid = spacegrids[0];

  double summed_grid = 0;
  // grid memory layout is (W)eight (T) Kinetic (V) potential
  for (int i = 0; i < 16000; i++)
    summed_grid += *(grid.getDataVector().begin() + i * 3 + 2) + *(grid.getDataVector().begin() + i * 3 + 1);

  // We can't just do this because you do not reduce values down to the head rank per particle logger
  // auto& pph_logger = dynamic_cast<PerParticleHamiltonianLogger&>(operator_ests[1].get());
  // so we just get this ranks;
  auto& pph_logger  = dynamic_cast<PerParticleHamiltonianLogger&>(crowd_operator_ests[0][1].get());
  auto expected_sum = pph_logger.sumOverAll();
  //Here we check the sum of logged energies against the total energy in the grid.

  // The head rank sees the reduction across ranks of the space grid. so its expectation is * num_ranks
  if (comm->rank() == 0)
    expected_sum *= comm->size();

  CHECK(summed_grid == Approx(expected_sum));
}

} // namespace qmcplusplus
