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


#ifndef QMCPLUSPLUS_ESTIMATORMANAGERBASETEST_HPP
#define QMCPLUSPLUS_ESTIMATORMANAGERBASETEST_HPP

#include "Estimators/EstimatorManagerBase.h"
#include "Estimators/tests/FakeEstimator.h"

class Communicate;

namespace qmcplusplus
{
namespace testing
{

/** Testing class breaking EstimatorManagerBase encapsultation
 *
 *  Wraps EstimatorManagerBase
 */
class EstimatorManagerBaseTest
{
public:
  EstimatorManagerBaseTest(Communicate* comm, int ranks);
  /** Quickly add scalar samples using FakeEstimator mock estimator. */
  void fakeSomeScalarSamples();
  /** call private EMB method and colelct EMBTs estimators_ */
  void collectScalarEstimators();
  /** for mpi test (it's trivial for 1 rank)
   *
   * only used by test_manager_mpi.cpp so implemented there.  
   */
  bool testMakeBlockAverages();
  EstimatorManagerBase em;
private:
  Communicate* comm_;
  std::vector<FakeEstimator> estimators_;

};

}
}

#endif /* QMCPLUSPLUS_ESTIMATORMANAGERBASETEST_HPP */
