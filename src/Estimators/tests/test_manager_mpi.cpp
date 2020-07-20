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

#include "catch.hpp"
#include "Message/Communicate.h"

#include "Platforms/Host/OutputManager.h"

#include "Estimators/EstimatorManagerBase.h"
#include "Estimators/tests/EstimatorManagerBaseTest.h"

namespace qmcplusplus
{

namespace testing {

bool EstimatorManagerBaseTest::testMakeBlockAverages()
{
  if(em.myComm->rank() == 1) {
    estimators_[1].scalars[0](3.0);
    estimators_[1].scalars[1](3.0);
    estimators_[1].scalars[2](3.0);
    estimators_[1].scalars[3](3.0);
  }

  // manipulation of state to arrive at to be tested state.
  // - From EstimatorManagerBase::reset
  em.weightInd = em.BlockProperties.add("BlockWeight");
  em.cpuInd    = em.BlockProperties.add("BlockCPU");
  em.acceptInd = em.BlockProperties.add("AcceptRatio");

  // - From EstimatorManagerBase::start
  em.PropertyCache.resize(em.BlockProperties.size());

  // - From EstimatorManagerBase::stopBlocknew
  //   three estimators
  //   - 2 with 1 sample 1
  //   - 1 with 2
  double block_weight = 0;
  std::for_each(estimators_.begin(),estimators_.end(),[&block_weight](auto& est){
                                                       block_weight += est.scalars[0].count();
                                                     });
  em.PropertyCache[em.weightInd] = block_weight;
  em.PropertyCache[em.cpuInd]    = 1.0;
  em.PropertyCache[em.acceptInd] = 1.0;

  RefVector<ScalarEstimatorBase> est_list = makeRefVector<ScalarEstimatorBase>(estimators_);
  em.collectScalarEstimators(est_list);

  em.makeBlockAverages();
  return true;
}

}

TEST_CASE("EstimatorManagerBase::makeBlockAverages()", "[estimators]")
{
  Communicate* c = OHMMS::Controller;
  int num_ranks = c->size();
  testing::EstimatorManagerBaseTest embt(c, num_ranks);

  embt.fakeSomeScalarSamples();
  embt.testMakeBlockAverages();

  // right now only rank() == 0 gets the actual averages
  if(c->rank() == 0)
  {
    double correct_value = ( 5.0 * num_ranks + 3.0) / (4 * (num_ranks -1) + 5);
    CHECK(embt.em.get_AverageCache()[0] == Approx(correct_value));
    correct_value = (8.0 * num_ranks + 3.0 ) / (4 * (num_ranks -1) + 5);
    CHECK(embt.em.get_AverageCache()[1] == Approx(correct_value));
    correct_value = (11.0 * num_ranks + 3.0 ) / (4 * (num_ranks -1) + 5);
    CHECK(embt.em.get_AverageCache()[2] == Approx(correct_value));
    correct_value = (14.0 * num_ranks + 3.0 ) / (4 * (num_ranks -1) + 5);
    CHECK(embt.em.get_AverageCache()[3] == Approx(correct_value));
  }
}

}
