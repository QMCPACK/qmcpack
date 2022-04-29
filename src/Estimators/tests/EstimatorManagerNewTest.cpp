//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "EstimatorManagerNewTest.h"
#include "Estimators/ScalarEstimatorBase.h"
#include "Platforms/Host/OutputManager.h"
#include "Estimators/tests/FakeEstimator.h"
#include "FakeOperatorEstimator.h"

namespace qmcplusplus
{
namespace testing
{
EstimatorManagerNewTest::EstimatorManagerNewTest(Communicate* comm, int ranks) : em(comm), comm_(comm)
{
  int num_ranks = comm_->size();
  if (num_ranks != ranks)
    throw std::runtime_error("Bad Rank Count, test expects different number of ranks.");

  app_log() << "running on " << num_ranks << '\n';
}

bool EstimatorManagerNewTest::testAddGetEstimator()
{
  // Must create on heap since the EstimatorManager destructor deletes all estimators
  auto fake_est_uptr      = std::make_unique<FakeEstimator>();
  FakeEstimator* fake_est = fake_est_uptr.get();

  em.add(std::move(fake_est_uptr), "fake");

  ScalarEstimatorBase* est2 = em.getEstimator("fake");
  FakeEstimator* fake_est2  = dynamic_cast<FakeEstimator*>(est2);

  return fake_est2 == fake_est;
}

void EstimatorManagerNewTest::fakeSomeScalarSamples()
{
  FakeEstimator fake_estimator;
  fake_estimator.scalars.resize(4);
  fake_estimator.scalars_saved.resize(4);
  fake_estimator.scalars[0](1.0);
  fake_estimator.scalars[1](2.0);
  fake_estimator.scalars[2](3.0);
  fake_estimator.scalars[3](4.0);

  //three estimators
  estimators_.push_back(fake_estimator);
  estimators_.push_back(fake_estimator);
  estimators_.push_back(fake_estimator);

  estimators_[2].scalars[0](2.0);
  estimators_[2].scalars[1](2.0);
  estimators_[2].scalars[2](2.0);
  estimators_[2].scalars[3](2.0);

  em.get_AverageCache().resize(4);
}

void EstimatorManagerNewTest::fakeSomeOperatorEstimatorSamples(int rank)
{
  em.operator_ests_.emplace_back(new FakeOperatorEstimator(comm_->size(), DataLocality::crowd));
  FakeOperatorEstimator& foe        = dynamic_cast<FakeOperatorEstimator&>(*(em.operator_ests_.back()));
  std::vector<QMCT::RealType>& data = foe.get_data();
  for (int id = 0; id < data.size(); ++id)
  {
    if (id > rank)
      data[id] += rank + 1;
  }
  foe.set_walker_weights(1);
}

std::vector<QMCTraits::RealType> EstimatorManagerNewTest::generateGoodOperatorData(int num_ranks)
{
  std::vector<QMCT::RealType> data(num_ranks * 10, 0.0);
  for (int ir = 0; ir < num_ranks; ++ir)
  {
    for (int id = 0; id < data.size(); ++id)
    {
      if (id > ir)
        data[id] += ir + 1;
    }
  }
  return data;
}

void EstimatorManagerNewTest::collectScalarEstimators()
{
  RefVector<ScalarEstimatorBase> est_list = makeRefVector<ScalarEstimatorBase>(estimators_);
  em.collectScalarEstimators(est_list);
}

void EstimatorManagerNewTest::testReduceOperatorEstimators() { em.reduceOperatorEstimators(); }

} // namespace testing
} // namespace qmcplusplus
