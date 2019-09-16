//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewin@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Estimators/EstimatorManagerBase.h"
#include "Estimators/tests/FakeEstimator.h"
#include "Estimators/ScalarEstimatorBase.h"


#include <stdio.h>
#include <sstream>

namespace qmcplusplus
{
TEST_CASE("EstimatorManagerBase", "[estimators]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate* c = OHMMS::Controller;

  EstimatorManagerBase em(c);

  REQUIRE(em.size() == 0);

  // Must create on heap since the EstimatorManager destructor deletes all estimators
  FakeEstimator* fake_est = new FakeEstimator;

  em.add(fake_est, "fake");

  ScalarEstimatorBase* est2 = em.getEstimator("fake");
  FakeEstimator* fake_est2  = dynamic_cast<FakeEstimator*>(est2);
  REQUIRE(fake_est2 != NULL);
  REQUIRE(fake_est2 == fake_est);

  // Check the copy constructor
  EstimatorManagerBase em2(em);
  REQUIRE(em.size() == 1);

  em.start(2, true);

  em.stop();
  // compute averages over threads
  //em.stop();
  em.reset();
}


TEST_CASE("EstimatorManagerBase::collectScalarEstimators", "[estimators]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate* c = OHMMS::Controller;

  FakeEstimator fake_estimator;
  fake_estimator.scalars.resize(4);
  fake_estimator.scalars_saved.resize(4);
  fake_estimator.scalars[0](1.0);
  fake_estimator.scalars[1](2.0);
  fake_estimator.scalars[2](3.0);
  fake_estimator.scalars[3](4.0);

  std::vector<FakeEstimator> estimators;
  estimators.push_back(fake_estimator);
  estimators.push_back(fake_estimator);
  estimators.push_back(fake_estimator);

  estimators[2].scalars[0](2.0);
  estimators[2].scalars[1](2.0);
  estimators[2].scalars[2](2.0);
  estimators[2].scalars[3](2.0);

  RefVector<ScalarEstimatorBase> est_list = makeRefVector<ScalarEstimatorBase>(estimators);
  EstimatorManagerBase em(c);

  em.get_AverageCache().resize(4);
  em.get_SquaredAverageCache().resize(4);
  est_list[0].get().scalars.resize(4);


  em.collectScalarEstimators(est_list);

  // We have three "walkers" two have one sample each of 1.0
  // one has two samples of 1.0, 2,0
  // Expected behavior is samples in a block for a walker are averaged
  // Walkers for blocks are averaged.  Walkers have equal weight.
  double correct_value = (2.0 + ((1.0 + 2.0) / 2)) / 3;
  REQUIRE(em.get_AverageCache()[0] == Approx(correct_value));
}




TEST_CASE("Estimator adhoc addVector operator", "[estimators]")
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
