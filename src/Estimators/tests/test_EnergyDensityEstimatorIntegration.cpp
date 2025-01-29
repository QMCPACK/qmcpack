//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "test_EnergyDensityEstimatorIntegration.h"
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
#include <for_testing/NativeInitializerPrint.hpp>
#include "GenerateRandomParticleSets.h"
#include "EstimatorManagerNewTest.h"
#include "EDenEstimatorManagerIntegrationTest.h"


constexpr bool generate_test_data = false;

namespace qmcplusplus
{
using MCPWalker = Walker<QMCTraits, PtclOnLatticeTraits>;

TEST_CASE("EnergyDensityEstimatorIntegration::multirank_reduction", "[estimators]")
{
  Communicate* comm = OHMMS::Controller;

  int num_walkers = 4;
  testing::EnergyDensityTest eden_test(comm, num_walkers, generate_test_data);
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

  emn.startDriverRun();

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

  auto operator_ests                  = emnta.getOperatorEstimators();
  NEEnergyDensityEstimator& e_den_est = dynamic_cast<NEEnergyDensityEstimator&>(operator_ests[0].get());

  auto spacegrids = e_den_est.getSpaceGrids();

  decltype(spacegrids)::value_type::type& grid = spacegrids[0];

  double summed_grid = 0;
  // grid memory layout is (W)eight (T) Kinetic (V) potential
  for (int i = 0; i < grid.nDomains(); i++)
    summed_grid += *(grid.getDataVector().begin() + i * 3 + 2) + *(grid.getDataVector().begin() + i * 3 + 1);

  auto block_weight = emc.get_block_weight();
  auto accept       = num_walkers * pset_list[0].getTotalNum();

  // We can't just do this because you do not reduce values down to the head rank per particle logger
  // auto& pph_logger = dynamic_cast<PerParticleHamiltonianLogger&>(operator_ests[1].get());
  // so we just get this ranks;
  auto& pph_logger = dynamic_cast<PerParticleHamiltonianLogger&>(crowd_operator_ests[0][1].get());
  using namespace std::string_literals;
  auto expected_sum = pph_logger.sumOverSome({"local_potential"s, "kinetic_energy"s, "ion_potential"s});
  //Here we check the sum of logged energies against the total energy in the grid.

  std::cout << "summed_grid: " << summed_grid << "  expected_sum: " << expected_sum << '\n';

  CHECK(summed_grid == Approx(expected_sum));
}

TEST_CASE("EnergyDensityEstimatorIntegration::operator_reporting", "[estimators]")
{
  Communicate* comm = OHMMS::Controller;
  int num_walkers   = 4;
  testing::EnergyDensityTest eden_test(comm, num_walkers, &testing::makeGoldWalkerElementsWithEEEIPS,
                                       generate_test_data);
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

  emn.startDriverRun();
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

  auto savePropertiesIntoWalker = [](QMCHamiltonian& ham, MCPWalker& walker) {
    ham.saveProperty(walker.getPropertyBase());
  };
  for (int iw = 0; iw < num_walkers; ++iw)
    savePropertiesIntoWalker(ham_list[iw], walker_list[iw]);

#ifndef NDEBUG
  {
    std::vector<RefVector<OperatorEstBase>> crowd_operator_ests = {emc.get_operator_estimators()};
    auto& e_den_est_crowd = dynamic_cast<NEEnergyDensityEstimator&>((crowd_operator_ests[0])[0].get());
    e_den_est_crowd.openDebugFile({"eden_values_" + std::to_string(comm->rank()) + "_particle.dat"});
  }
#endif

  emc.accumulate(walker_list, pset_list, twf_list, ham_list, rng);

  std::vector<RefVector<OperatorEstBase>> crowd_operator_ests = {emc.get_operator_estimators()};
  emn.collectOperatorEstimators(crowd_operator_ests);

  RefVector<ScalarEstimatorBase> main_scalar_estimators;
  main_scalar_estimators.push_back(emc.get_main_estimator());
  emn.collectMainEstimators(main_scalar_estimators);

  testing::EstimatorManagerNewTestAccess emnta(emn);
  auto operator_ests = emnta.getOperatorEstimators();

  NEEnergyDensityEstimator& e_den_est = dynamic_cast<NEEnergyDensityEstimator&>(operator_ests[0].get());

  auto spacegrids = e_den_est.getSpaceGrids();

  decltype(spacegrids)::value_type::type& grid = spacegrids[0];

  double summed_grid = 0;
  // grid memory layout is (W)eight (T) Kinetic (V) potential
  for (int i = 0; i < grid.nDomains(); i++)
  {
    summed_grid += *(grid.getDataVector().begin() + i * 3 + 2) + *(grid.getDataVector().begin() + i * 3 + 1);
  }
  // We can't just do this because you do not reduce values down to the head rank per particle logger
  // auto& pph_logger = dynamic_cast<PerParticleHamiltonianLogger&>(operator_ests[1].get());
  // so we just get this ranks;
  auto& pph_logger = dynamic_cast<PerParticleHamiltonianLogger&>(crowd_operator_ests[0][1].get());
  using namespace std::string_literals;
  auto expected_sum = pph_logger.sumOverSome({"local_potential"s, "kinetic_energy"s, "ion_potential"s});
  //Here we check the sum of logged energies against the total energy in the grid.

  auto block_weight = emc.get_block_weight();
  auto accept       = num_walkers * pset_list[0].getTotalNum();

  CHECK(summed_grid == Approx(expected_sum));
  auto debug_sum = testing::cannedSum();
  CHECK(summed_grid == Approx(debug_sum));
}


TEST_CASE("EnergyDensityEstimatorIntegration::normalization", "[estimators]")
{
  int num_walkers   = 4;
  Communicate* comm = OHMMS::Controller;
  testing::EDenEstimatorManagerIntegrationTest eden_emn_integration_test(comm, num_walkers);
  auto walker_list = eden_emn_integration_test.getWalkerList();
  auto ham_list    = eden_emn_integration_test.getHamList();

  auto savePropertiesIntoWalker = [](QMCHamiltonian& ham, MCPWalker& walker) {
    ham.saveProperty(walker.getPropertyBase());
  };
  for (int iw = 0; iw < num_walkers; ++iw)
    savePropertiesIntoWalker(ham_list[iw], walker_list[iw]);

  auto& emc = eden_emn_integration_test.getEmc();

  std::vector<RefVector<OperatorEstBase>> crowd_operator_ests = {emc.get_operator_estimators()};

  // We can't just do this because you do not reduce values down to the head rank per particle logger
  // auto& pph_logger = dynamic_cast<PerParticleHamiltonianLogger&>(operator_ests[1].get());
  // so we just get this ranks;
  auto& pph_logger = dynamic_cast<PerParticleHamiltonianLogger&>(crowd_operator_ests[0][1].get());
  using namespace std::string_literals;
  auto expected_sum = pph_logger.sumOverSome({"local_potential"s, "kinetic_energy"s, "ion_potential"s});

  std::cout << "expected_sum: " << expected_sum << '\n';

  auto pset_list = eden_emn_integration_test.getPSetList();
  auto twf_list  = eden_emn_integration_test.getTwfList();
  auto& rng      = eden_emn_integration_test.getRng();
  // The logger writes out on accumualtes so the sum must be acquired before
  emc.accumulate(walker_list, pset_list, twf_list, ham_list, rng);

  NEEnergyDensityEstimator& e_den_est = dynamic_cast<NEEnergyDensityEstimator&>(crowd_operator_ests[0][0].get());
  auto spacegrids                     = e_den_est.getSpaceGrids();
  decltype(spacegrids)::value_type::type& grid = spacegrids[0];
  double summed_grid                           = 0;
  // grid memory layout is (W)eight (T) Kinetic (V) potential
  for (int i = 0; i < grid.nDomains(); i++)
  {
    summed_grid += *(grid.getDataVector().begin() + i * 3 + 2) + *(grid.getDataVector().begin() + i * 3 + 1);
  }

  CHECK(summed_grid == Approx(expected_sum));
  auto debug_sum = testing::cannedSum();
  CHECK(summed_grid == Approx(debug_sum));


  std::cout << "summed grid: " << summed_grid << '\n';

  int num_steps = 4;
  
  double expected_sum2 = 0.0;
  for (int step = 1; step < num_steps; ++step)
  {
    // later steps
    auto& gold_elements = eden_emn_integration_test.getGoldElements();
    testing::particleSetsToRandomPositions(pset_list, gold_elements.pset_ions,
                                           eden_emn_integration_test.deterministic_rs_2, generate_test_data);

    eden_emn_integration_test.updateAndEvaluate();

    expected_sum2 += pph_logger.sumOverSome({"local_potential"s, "kinetic_energy"s, "ion_potential"s});

    std::cout << "expected sum2: " << expected_sum2 << '\n';
    emc.accumulate(walker_list, pset_list, twf_list, ham_list, rng);

    summed_grid = 0;
    for (int i = 0; i < grid.nDomains(); i++)
    {
      summed_grid += *(grid.getDataVector().begin() + i * 3 + 2) + *(grid.getDataVector().begin() + i * 3 + 1);
    }

    CHECK(summed_grid == Approx(expected_sum + expected_sum2));
  }
  
  auto block_weight = emc.get_block_weight();
  std::cout << "block weight: " << block_weight << '\n';
  auto accept = num_walkers * pset_list[0].getTotalNum();

  RefVector<ScalarEstimatorBase> main_scalar_estimators;
  main_scalar_estimators.push_back(emc.get_main_estimator());
  auto& emn = eden_emn_integration_test.getEmn();
  emn.collectMainEstimators(main_scalar_estimators);
  emn.collectOperatorEstimators(crowd_operator_ests);
  testing::EstimatorManagerNewTestAccess emnta(emn);

  auto operator_ests                               = emnta.getOperatorEstimators();
  NEEnergyDensityEstimator& app_e_den_est          = dynamic_cast<NEEnergyDensityEstimator&>(operator_ests[0].get());
  auto app_spacegrids                              = app_e_den_est.getSpaceGrids();
  decltype(spacegrids)::value_type::type& app_grid = app_spacegrids[0];

  summed_grid = 0;
  // grid memory layout is (W)eight (T) Kinetic (V) potential
  for (int i = 0; i < app_grid.nDomains(); i++)
  {
    summed_grid += *(app_grid.getDataVector().begin() + i * 3 + 2) + *(app_grid.getDataVector().begin() + i * 3 + 1);
  }

  std::cout << "summed grid: " << summed_grid << '\n';
  CHECK(summed_grid == Approx(expected_sum + expected_sum2));
  std::cout << "walkers weight: " << app_e_den_est.get_walkers_weight() << '\n';

  emnta.stopBlockUpToWrite(accept, 0, block_weight);
  std::cout << "walkers weight after reduction: " << app_e_den_est.get_walkers_weight() << '\n';

  if (comm->rank() == 0)
  {
    summed_grid = 0;
    // grid memory layout is (W)eight (T) Kinetic (V) potential

    for (int i = 0; i < app_grid.nDomains(); i++)
    {
      summed_grid += *(app_grid.getDataVector().begin() + i * 3 + 2) + *(app_grid.getDataVector().begin() + i * 3 + 1);
    }
    std::cout << "summed grid: " << summed_grid << '\n';

    // This is by 8 because expected_sum + expected_sum2 are only the local grids
    CHECK(summed_grid == Approx((expected_sum + expected_sum2) / (num_walkers * num_steps)));
  }
}

} // namespace qmcplusplus
