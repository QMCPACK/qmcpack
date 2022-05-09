//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "Message/Communicate.h"

#include "Platforms/Host/OutputManager.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Estimators/EstimatorManagerNew.h"
#include "Estimators/tests/EstimatorManagerNewTest.h"

namespace qmcplusplus
{
namespace testing
{
bool EstimatorManagerNewTest::testMakeBlockAverages()
{
  if (em.my_comm_->rank() == 1)
  {
    estimators_[1].scalars[0](3.0);
    estimators_[1].scalars[1](3.0);
    estimators_[1].scalars[2](3.0);
    estimators_[1].scalars[3](3.0);
  }

  // manipulation of state to arrive at to be tested state.
  // - From EstimatorManagerBase::reset
  em.weightInd      = em.BlockProperties.add("BlockWeight");
  em.cpuInd         = em.BlockProperties.add("BlockCPU");
  em.acceptRatioInd = em.BlockProperties.add("AcceptRatio");

  // - From EstimatorManagerBase::start
  em.PropertyCache.resize(em.BlockProperties.size());

  // - From EstimatorManagerBase::stopBlocknew
  //   three estimators
  //   - 2 with 1 sample 1
  //   - 1 with 2
  double block_weight = 0;
  std::for_each(estimators_.begin(), estimators_.end(),
                [&block_weight](auto& est) { block_weight += est.scalars[0].count(); });
  em.PropertyCache[em.weightInd] = block_weight;
  em.PropertyCache[em.cpuInd]    = 1.0;

  RefVector<ScalarEstimatorBase> est_lists(makeRefVector<ScalarEstimatorBase>(estimators_));
  em.collectMainEstimators(est_lists);

  unsigned long accepts = 4;
  unsigned long rejects = 1;
  em.makeBlockAverages(accepts, rejects);
  return true;
}

} // namespace testing

TEST_CASE("EstimatorManagerNew::makeBlockAverages()", "[estimators]")
{
  Communicate* c = OHMMS::Controller;
  int num_ranks  = c->size();
  QMCHamiltonian ham;
  testing::EstimatorManagerNewTest embt(ham, c, num_ranks);

  embt.fakeMainScalarSamples();
  embt.testMakeBlockAverages();

  // right now only rank() == 0 gets the actual averages
  if (c->rank() == 0)
  {
    double correct_value = (5.0 * num_ranks + 3.0) / (4 * (num_ranks - 1) + 5);
    CHECK(embt.em.get_AverageCache()[0] == Approx(correct_value));
    correct_value = (8.0 * num_ranks + 3.0) / (4 * (num_ranks - 1) + 5);
    CHECK(embt.em.get_AverageCache()[1] == Approx(correct_value));
    correct_value = (11.0 * num_ranks + 3.0) / (4 * (num_ranks - 1) + 5);
    CHECK(embt.em.get_AverageCache()[2] == Approx(correct_value));
    correct_value = (14.0 * num_ranks + 3.0) / (4 * (num_ranks - 1) + 5);
    CHECK(embt.em.get_AverageCache()[3] == Approx(correct_value));
  }
}

TEST_CASE("EstimatorManagerNew::reduceOperatorestimators()", "[estimators]")
{
  Communicate* c = OHMMS::Controller;
  int num_ranks  = c->size();
  QMCHamiltonian ham;
  testing::EstimatorManagerNewTest embt(ham, c, num_ranks);

  embt.fakeSomeOperatorEstimatorSamples(c->rank());
  std::vector<QMCTraits::RealType> good_data = embt.generateGoodOperatorData(num_ranks);

  // Normalization is done by reduceOperatorEstimators based on the the total weight of the
  // estimators for that block.
  embt.testReduceOperatorEstimators();

  if (c->rank() == 0)
  {
    auto& test_data = embt.get_operator_data();

    QMCTraits::RealType norm = 1.0 / static_cast<QMCTraits::RealType>(num_ranks);
    for (size_t i = 0; i < test_data.size(); ++i)
    {
      QMCTraits::RealType norm_good_data = good_data[i] * norm;
      if (norm_good_data != test_data[i])
      {
        FAIL_CHECK("norm_good_data " << norm_good_data << " != test_data " << test_data[i] << " at index " << i);
        break;
      }
    }
  }
}

} // namespace qmcplusplus
