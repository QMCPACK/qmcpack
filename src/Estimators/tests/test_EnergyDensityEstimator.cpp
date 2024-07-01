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

#include "EnergyDensityEstimator.h"
#include <iostream>
#include "EstimatorTesting.h"
#include "EnergyDensityTest.h"
#include "PerParticleHamiltonianLogger.h"
#include "ValidEnergyDensityInput.h"
#include "Particle/Walker.h"
#include "Particle/tests/MinimalParticlePool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
#include "QMCHamiltonians/tests/MinimalHamiltonianPool.h"
#include "Utilities/for_testing/NativeInitializerPrint.hpp"

constexpr bool generate_test_data = false;

namespace qmcplusplus
{

TEST_CASE("NEEnergyDensityEstimator::Constructor", "[estimators]")
{
  Communicate* comm;
  comm = OHMMS::Controller;

  ParticleSetPool particle_pool{MinimalParticlePool::make_diamondC_1x1x1(comm)};

  ParticleSet pset_elec{*(particle_pool.getParticleSet("e"))};
  ParticleSet pset_ions{*(particle_pool.getParticleSet("ion"))};

  pset_elec.R =
      ParticleSet::ParticlePos{{1.451870349, 1.381521229, 1.165202269}, {1.244515371, 1.382273176, 1.21105285},
                               {0.000459944, 1.329603408, 1.265030556}, {0.748660329, 1.63420622, 1.393637791},
                               {0.033228526, 1.391869137, 0.654413566}, {1.114198787, 1.654334594, 0.231075822},
                               {1.657151589, 0.883870516, 1.201243939}, {0.97317591, 1.245644974, 0.284564732}};

  Libxml2Document doc;
  using Input = testing::EnergyDensityInputs;
  Input input;
  bool okay       = doc.parseFromString(input[Input::valid::CELL]);
  xmlNodePtr node = doc.getRoot();
  EnergyDensityInput edein{node};
  {
    NEEnergyDensityEstimator e_den_est(edein, particle_pool.getPool());
  }
}

TEST_CASE("NEEnergyDensityEstimator::spawnCrowdClone", "[estimators]")
{
  Communicate* comm;
  comm = OHMMS::Controller;

  ParticleSetPool particle_pool{MinimalParticlePool::make_diamondC_1x1x1(comm)};

  ParticleSet pset_elec{*(particle_pool.getParticleSet("e"))};
  ParticleSet pset_ions{*(particle_pool.getParticleSet("ion"))};

  pset_elec.R =
      ParticleSet::ParticlePos{{1.451870349, 1.381521229, 1.165202269}, {1.244515371, 1.382273176, 1.21105285},
                               {0.000459944, 1.329603408, 1.265030556}, {0.748660329, 1.63420622, 1.393637791},
                               {0.033228526, 1.391869137, 0.654413566}, {1.114198787, 1.654334594, 0.231075822},
                               {1.657151589, 0.883870516, 1.201243939}, {0.97317591, 1.245644974, 0.284564732}};

  Libxml2Document doc;
  using Input = testing::EnergyDensityInputs;
  Input input;
  bool okay       = doc.parseFromString(input[Input::valid::CELL]);
  xmlNodePtr node = doc.getRoot();
  EnergyDensityInput edein{node};
  {
    NEEnergyDensityEstimator original_e_den_est(edein, particle_pool.getPool());

    auto clone = original_e_den_est.spawnCrowdClone();
    REQUIRE(clone != nullptr);
    REQUIRE(clone.get() != &original_e_den_est);
    REQUIRE(dynamic_cast<decltype(&original_e_den_est)>(clone.get()) != nullptr);
  }
}

TEST_CASE("NEEnergyDensityEstimator::AccumulateIntegration", "[estimators]")
{
  Communicate* comm = OHMMS::Controller;

  testing::EnergyDensityTest eden_test(comm, 4 /*num_walkers*/, generate_test_data);
  
  // using MCPWalker = typename OperatorEstBase::MCPWalker;
  // std::vector<OperatorEstBase::MCPWalker> walkers(num_walkers, MCPWalker(gold_elem.pset_elec.getTotalNum()));

  // auto walker_refs = makeRefVector<MCPWalker>(walkers);
  // RefVectorWithLeader<MCPWalker> walker_list{walker_refs[0], walker_refs};

  // RefVector<QMCHamiltonian> ham_refs = convertUPtrToRefVector(hams);

  // RefVectorWithLeader<QMCHamiltonian> ham_list{ham_refs[0], ham_refs};
  // ResourceCollection ham_res("test_ham_res");
  // ham_list.getLeader().createResource(ham_res);
  // ResourceCollectionTeamLock<QMCHamiltonian> ham_lock(ham_res, ham_list);
  auto ham_list = eden_test.getHamList();
  auto& ham_leader = ham_list.getLeader();
  auto ham_lock = eden_test.makeTeamLock(eden_test.getHamRes(), ham_list);
  eden_test.getEnergyDensityEstimator().registerListeners(ham_leader);

  PerParticleHamiltonianLogger pph_logger({}, 0);
  pph_logger.registerListeners(ham_leader);

  auto pset_list = eden_test.getPSetList();
  auto pset_lock = eden_test.makeTeamLock(eden_test.getPSetRes(), pset_list);

  pset_list[0].L[0] = 1.0;
  pset_list[1].L[1] = 1.0;
  pset_list[2].L[2] = 1.0;
  pset_list[3].L[3] = 1.0;

  ParticleSet::mw_update(pset_list);
  auto twf_list = eden_test.getTwfList();
  auto twf_lock = eden_test.makeTeamLock(eden_test.getTwfRes(), twf_list);
  TrialWaveFunction::mw_evaluateLog(twf_list, pset_list);
  QMCHamiltonian::mw_evaluate(ham_list, twf_list, pset_list);

  hdf_archive hd;
  std::string test_file{"ede_test.hdf"};
  bool okay = hd.create(test_file);
  REQUIRE(okay);
  std::vector<ObservableHelper> h5desc;
  auto& e_den_est = eden_test.getEnergyDensityEstimator();
  e_den_est.registerOperatorEstimator(hd);

  StdRandom<double> rng;
  rng.init(101);

  e_den_est.accumulate(eden_test.getWalkerList(), pset_list, twf_list, ham_list, rng);
  auto spacegrids = e_den_est.getSpaceGrids();

  decltype(spacegrids)::value_type::type& grid = spacegrids[0];

  double summed_grid = 0;
  // grid memory layout is (W)eight (T) Kinetic (V) potential
  for (int i = 0; i < 16000; i++)
    summed_grid += *(grid.getDataVector().begin() + i * 3 + 2) + *(grid.getDataVector().begin() + i * 3 + 1);

  auto expected_sum = pph_logger.sumOverAll();
  //Here we check the sum of logged energies against the total energy in the grid.
  CHECK(summed_grid == Approx(expected_sum));

  e_den_est.write(hd);
}

TEST_CASE("NEEnergyDensityEstimator::Collect", "[estimators]") {}

} // namespace qmcplusplus
