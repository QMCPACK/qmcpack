//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, , doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_QMCDRIVERNEWTESTWRAPPER_H
#define QMCPLUSPLUS_QMCDRIVERNEWTESTWRAPPER_H
#include "QMCDrivers/QMCDriverNew.h"
#include "QMCDrivers/DriverTraits.h"
#include "Particle/SampleStack.h"
#include "Concurrency/ParallelExecutor.hpp"
#include "Message/UniformCommunicateError.h"

namespace qmcplusplus
{
/** Unit testing an impure virtual base class
 *  requires a absolute minimal subtype
 */
namespace testing
{
class QMCDriverNewTestWrapper : public QMCDriverNew
{
public:
  using Base = QMCDriverNew;
  QMCDriverNewTestWrapper(QMCDriverInput&& input,
                          MCPopulation&& population,
                          TrialWaveFunction& psi,
                          QMCHamiltonian& h,
                          SampleStack samples,
                          Communicate* comm)
      : QMCDriverNew(std::move(input),
                     std::move(population),
                     psi,
                     h,
                     "QMCDriverTestWrapper::",
                     comm,
                     "QMCDriverNewTestWrapper")
  {}

  ~QMCDriverNewTestWrapper() {}

  QMCRunType getRunType() { return QMCRunType::DUMMY; }

  void process(xmlNodePtr node)
  {
    // We want to test the reserve ability as well
    AdjustedWalkerCounts awc =
        adjustGlobalWalkerCount(myComm->size(), myComm->rank(), qmcdriver_input_.get_total_walkers(),
                                qmcdriver_input_.get_walkers_per_rank(), 1.0, qmcdriver_input_.get_num_crowds());

    Base::startup(node, awc);
  }

  void testAdjustGlobalWalkerCount()
  {
    QMCDriverNew::AdjustedWalkerCounts awc = adjustGlobalWalkerCount(4, 1, 64, 0, 1.0, 0);
    CHECK(awc.global_walkers == 64);
    CHECK(awc.walkers_per_crowd.size() == 8);
    CHECK(awc.walkers_per_rank[0] == 16);
    CHECK(awc.walkers_per_rank[3] == 16);
    CHECK(awc.walkers_per_crowd[4] == 2);
    CHECK(awc.walkers_per_crowd[7] == 2);

    awc = adjustGlobalWalkerCount(4, 1, 63, 0, 1.0, 4);
    CHECK(awc.global_walkers == 63);
    CHECK(awc.walkers_per_crowd.size() == 4);
    CHECK(awc.walkers_per_rank[0] == 16);
    CHECK(awc.walkers_per_rank[3] == 15);
    CHECK(awc.walkers_per_crowd[0] == 4);
    CHECK(awc.walkers_per_crowd[3] == 4);

    awc = adjustGlobalWalkerCount(4, 3, 63, 0, 1.0, 4);
    CHECK(awc.global_walkers == 63);
    CHECK(awc.walkers_per_crowd.size() == 4);
    CHECK(awc.walkers_per_rank[0] == 16);
    CHECK(awc.walkers_per_rank[3] == 15);
    CHECK(awc.walkers_per_crowd[0] == 4);
    CHECK(awc.walkers_per_crowd[3] == 3);

    awc = adjustGlobalWalkerCount(4, 3, 0, 32, 1.0, 4);
    CHECK(awc.global_walkers == 128);
    CHECK(awc.walkers_per_crowd.size() == 4);
    CHECK(awc.walkers_per_rank[0] == 32);
    CHECK(awc.walkers_per_rank[3] == 32);
    CHECK(awc.walkers_per_crowd[0] == 8);
    CHECK(awc.walkers_per_crowd[3] == 8);

    awc = adjustGlobalWalkerCount(1, 0, 7, 0, 1.0, 8);
    CHECK(awc.global_walkers == 7);
    CHECK(awc.walkers_per_crowd.size() == 8);
    CHECK(awc.walkers_per_rank.size() == 1);
    CHECK(awc.walkers_per_rank[0] == 7);
    CHECK(awc.walkers_per_crowd[0] == 1);
    CHECK(awc.walkers_per_crowd[7] == 0);

    awc = adjustGlobalWalkerCount(2, 0, 28, 0, 1.0, 0);
    CHECK(awc.global_walkers == 28);
    CHECK(awc.walkers_per_crowd.size() == Concurrency::maxCapacity());
    CHECK(awc.walkers_per_rank.size() == 2);
    CHECK(awc.walkers_per_rank[0] == 14);
    // \todo for std::thread these will be ones
    CHECK(awc.walkers_per_crowd[0] == 2);
    CHECK(awc.walkers_per_crowd[7] == 1);

    awc = adjustGlobalWalkerCount(2, 1, 27, 0, 1.0, 4);
    CHECK(awc.global_walkers == 27);
    CHECK(awc.walkers_per_crowd.size() == 4);
    CHECK(awc.walkers_per_rank.size() == 2);
    CHECK(awc.walkers_per_rank[0] == 14);
    CHECK(awc.walkers_per_rank[1] == 13);
    CHECK(awc.walkers_per_crowd[0] == 4);
    CHECK(awc.walkers_per_crowd[3] == 3);

    // Ask for 27 total walkers on 2 ranks of 11 walkers (inconsistent input)
    // results in fatal exception on all ranks.
    CHECK_THROWS_AS(adjustGlobalWalkerCount(2, 1, 27, 11, 1.0, 4), UniformCommunicateError);
    // Ask for 14 total walkers on 16 ranks (inconsistent input)
    // results in fatal exception on all ranks.
    CHECK_THROWS_AS(adjustGlobalWalkerCount(16, 0, 14, 0, 0, 0), UniformCommunicateError);
  }

  bool run() { return false; }

  int get_num_crowds() { return crowds_.size(); }

  template<class CONCURRENCY>
  struct TestNumCrowdsVsNumThreads
  {
    void operator()(int num_crowds);
  };

private:
};

template<class CONCURRENCY>
void QMCDriverNewTestWrapper::TestNumCrowdsVsNumThreads<CONCURRENCY>::operator()(int num_crowds)
{}

template<>
void QMCDriverNewTestWrapper::TestNumCrowdsVsNumThreads<ParallelExecutor<Executor::OPENMP>>::operator()(int num_crowds)
{
  if (Concurrency::maxCapacity<>() != 8)
    throw std::runtime_error("OMP_NUM_THREADS must be 8 for this test.");
  if (num_crowds > 8)
    CHECK_THROWS_AS(checkNumCrowdsLTNumThreads(num_crowds), UniformCommunicateError);
  else
    checkNumCrowdsLTNumThreads(num_crowds);
  return;
}


} // namespace testing
} // namespace qmcplusplus
#endif
