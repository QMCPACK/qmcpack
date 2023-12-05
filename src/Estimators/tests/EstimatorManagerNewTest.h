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


#ifndef QMCPLUSPLUS_ESTIMATORMANAGERNEWTEST_HPP
#define QMCPLUSPLUS_ESTIMATORMANAGERNEWTEST_HPP

#include "Estimators/EstimatorManagerNew.h"
#include "Estimators/tests/FakeEstimator.h"

class Communicate;

namespace qmcplusplus
{
namespace testing
{

/** Testing class breaking EstimatorManagerNew encapsultation
 *
 *  Wraps EstimatorManagerNew
 */
class EstimatorManagerNewTest
{
public:
  using QMCT = QMCTraits;
  
  EstimatorManagerNewTest(const QMCHamiltonian& ham, Communicate* comm, int ranks);
  /** Quickly add main scalar samples using FakeEstimator mock estimator. */
  void fakeMainScalarSamples();
  /** Quickly add scalar samples using FakeEstimator mock estimator. */
  void fakeScalarSamplesAndCollect();
  /** Quickly add scalar samples using FakeOperatorEstimator mock estimator. */
  void fakeSomeOperatorEstimatorSamples(int rank);
  /** call private EMB method and collect EMBTs estimators_ as main_estimators*/
  void collectMainEstimators();
  /** call private EMB method and colelct EMBTs estimators_ */
  void collectScalarEstimators();
  /** reduce the OperatorEstimators onto the EstimatorManagerNew copy. */
  void collectOperatorEstimators();
  /** for mpi test (it's trivial for 1 rank)
   *
   * only used by test_manager_mpi.cpp so implemented there.  
   */
  std::vector<QMCT::RealType> generateGoodOperatorData(int num_ranks);
  /// test replacing the main estimator
  bool testReplaceMainEstimator();
  
  bool testMakeBlockAverages();
  void testReduceOperatorEstimators();

  std::vector<QMCT::RealType>& get_operator_data() { return em.operator_ests_[0]->get_data(); }
  
  EstimatorManagerNew em;
private:
  Communicate* comm_;
  std::vector<FakeEstimator> estimators_;
  std::vector<RefVector<ScalarEstimatorBase>> scalar_estimators_;
};

class EstimatorManagerNewTestAccess {
public:
  EstimatorManagerNewTestAccess(EstimatorManagerNew& emn) : emn_(emn) {}
  const ScalarEstimatorBase& getMainEstimator() { return *(emn_.main_estimator_.get()); }
private:
  EstimatorManagerNew& emn_;
};

}
}

#endif /* QMCPLUSPLUS_ESTIMATORMANAGERNEWTEST_HPP */
