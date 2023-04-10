//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File refactored from: Refactored from test_manager.cpp
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Estimators/EstimatorManagerNew.h"
#include "Estimators/ScalarEstimatorBase.h"
#include "Estimators/tests/EstimatorManagerNewTest.h"
#include "Estimators/EstimatorInputDelegates.h"
#include "EstimatorManagerInputTest.h"
#include "Particle/tests/MinimalParticlePool.h"
#include "QMCHamiltonians/tests/MinimalHamiltonianPool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
#include "Utilities/ProjectData.h"
#include <stdio.h>
#include <sstream>

namespace qmcplusplus
{

TEST_CASE("EstimatorManagerNew::EstimatorManager(comm)", "[estimators]")
{
  Communicate* c = OHMMS::Controller;
  using namespace testing;

  QMCHamiltonian ham;
  testing::EstimatorManagerNewTest embt(ham, c, 1);

  // Just checking here that we have an empty estimator manager in embt.

  CHECK(embt.em.getNumEstimators() == 0);
  CHECK(embt.em.getNumScalarEstimators() == 0);
}

TEST_CASE("EstimatorManagerNew::EstimatorManagerNew(EstimatorManagerInput,...)", "[estimators]")
{
  ProjectData test_project("test", ProjectData::DriverVersion::BATCH);
  Communicate* comm = OHMMS::Controller;

  using namespace testing;
  Libxml2Document estimators_doc = createEstimatorManagerNewVMCInputXML();
  EstimatorManagerInput emi(estimators_doc.getRoot());

  CHECK(emi.get_estimator_inputs().size() == 2);
  CHECK(emi.get_scalar_estimator_inputs().size() == 1);

  auto particle_pool = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool =
      MinimalWaveFunctionPool::make_diamondC_1x1x1(test_project.getRuntimeOptions(), comm, particle_pool);
  auto& pset            = *(particle_pool.getParticleSet("e"));
  auto hamiltonian_pool = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);
  auto& twf             = *(wavefunction_pool.getWaveFunction("wavefunction"));
  auto& ham             = *(hamiltonian_pool.getPrimary());
  EstimatorManagerNew emn(comm, std::move(emi), ham, pset, twf);

  CHECK(emn.getNumEstimators() == 2);
  // Because the only scalar estimator becomes the main estimator.
  CHECK(emn.getNumScalarEstimators() == 0);
  EstimatorManagerNewTestAccess emnta(emn);
  CHECK(emnta.getMainEstimator().getName() == "LocalEnergyEstimator");


  // Check behavior when multiple "main" estimators are in the input
  // Each should override the previous main input as this is the behavior of legacy.
  Libxml2Document estimators_doc2 = createEstimatorManagerNewInputXML();
  EstimatorManagerInput emi2(estimators_doc2.getRoot());

  CHECK(emi2.get_estimator_inputs().size() == 2);
  CHECK(emi2.get_scalar_estimator_inputs().size() == 5);

  EstimatorManagerNew emn2(comm, std::move(emi2), ham, pset, twf);

  CHECK(emn2.getNumEstimators() == 2);
  // Because the only scalar estimator becomes the main estimator.
  CHECK(emn2.getNumScalarEstimators() == 0);
  EstimatorManagerNewTestAccess emnta2(emn2);
  CHECK(emnta2.getMainEstimator().getName() == "RMCLocalEnergyEstimator");
}

TEST_CASE("EstimatorManagerNew::collectMainEstimators", "[estimators]")
{
  Communicate* c = OHMMS::Controller;

  QMCHamiltonian ham;
  testing::EstimatorManagerNewTest embt(ham, c, 1);
  // by design we have done no averaging here
  // the division by total weight happens only when a block is over and the
  // accumulated data has been reduced down.  So here there should just be simple sums.

  REQUIRE(embt.testReplaceMainEstimator());

  embt.fakeMainScalarSamples();
  embt.collectMainEstimators();
  double correct_value = 5.0;
  CHECK(embt.em.get_AverageCache()[0] == Approx(correct_value));
  correct_value = 8.0;
  CHECK(embt.em.get_AverageCache()[1] == Approx(correct_value));
  correct_value = 11.0;
  CHECK(embt.em.get_AverageCache()[2] == Approx(correct_value));
  correct_value = 14.0;
  CHECK(embt.em.get_AverageCache()[3] == Approx(correct_value));
}

TEST_CASE("EstimatorManagerNew::collectScalarEstimators", "[estimators]")
{
  Communicate* c = OHMMS::Controller;

  QMCHamiltonian ham;
  testing::EstimatorManagerNewTest embt(ham, c, 1);
  // by design we have done no averaging here
  // the division by total weight happens only when a block is over and the
  // accumulated data has been reduced down.  So here there should just be simple sums.

  embt.fakeScalarSamplesAndCollect();
  double correct_value = 11.0;
  CHECK(embt.em.get_AverageCache()[0] == Approx(correct_value));
  correct_value = 20.0;
  CHECK(embt.em.get_AverageCache()[1] == Approx(correct_value));
  correct_value = 29.0;
  CHECK(embt.em.get_AverageCache()[2] == Approx(correct_value));
  correct_value = 38.0;
  CHECK(embt.em.get_AverageCache()[3] == Approx(correct_value));
}

TEST_CASE("EstimatorManagerNew adhoc addVector operator", "[estimators]")
{
  int num_scalars = 3;
  std::vector<double> vec_a{1.0, 2.0, 3.0};
  std::vector<double> vec_b{2.0, 3.0, 4.0};
  std::vector<std::vector<double>> est{vec_a, vec_b};
  auto addVectors = [](const auto& vec_a, const auto& vec_b) {
    std::vector<ScalarEstimatorBase::RealType> result_vector(vec_a.size(), 0.0);
    for (int i = 0; i < vec_a.size(); ++i)
      result_vector[i] = vec_a[i] + vec_b[i];
    return result_vector;
  };

  std::vector<double> reduced_scalars(num_scalars);
  reduced_scalars = std::accumulate(est.begin(), est.end(), std::vector<double>(num_scalars, 0.0), addVectors);
  std::vector<double> correct{3.0, 5.0, 7.0};
  REQUIRE(reduced_scalars == correct);
}

} // namespace qmcplusplus
