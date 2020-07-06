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
                          MCPopulation& population,
                          TrialWaveFunction& psi,
                          QMCHamiltonian& h,
                          WaveFunctionPool& ppool,
                          SampleStack samples,
                          Communicate* comm)
      : QMCDriverNew(std::move(input), population, psi, h, ppool, "QMCDriverTestWrapper::", comm)
  {
    QMCType = "QMCTesting";
  }

  ~QMCDriverNewTestWrapper() {}
  
  QMCRunType getRunType() { return QMCRunType::DUMMY; }

  void process(xmlNodePtr node)
  {
    // We want to test the reserve ability as well
    AdjustedWalkerCounts awc = adjustGlobalWalkerCount(myComm->size(), myComm->rank(), qmcdriver_input_.get_total_walkers(), qmcdriver_input_.get_walkers_per_rank(), 1.0, qmcdriver_input_.get_num_crowds());

    Base::startup(node, awc);
  }

  void testAdjustGlobalWalkerCount()
  {
    QMCDriverNew::AdjustedWalkerCounts awc = adjustGlobalWalkerCount(4,1,64,0,1.0,0);
    CHECK(awc.global_walkers == 64);
    CHECK(awc.walkers_per_crowd.size() == 8);
    CHECK(awc.walkers_per_rank[0] == 16);
    CHECK(awc.walkers_per_rank[3] == 16);
    CHECK(awc.walkers_per_crowd[4] == 2);
    CHECK(awc.walkers_per_crowd[7] == 2);

     awc = adjustGlobalWalkerCount(4,1,63,0,1.0,4);
    CHECK(awc.global_walkers == 63);
    CHECK(awc.walkers_per_crowd.size() == 4);
    CHECK(awc.walkers_per_rank[0] == 16);
    CHECK(awc.walkers_per_rank[3] == 15);
    CHECK(awc.walkers_per_crowd[0] == 4);
    CHECK(awc.walkers_per_crowd[3] == 4);

     awc = adjustGlobalWalkerCount(4,3,63,0,1.0,4);
    CHECK(awc.global_walkers == 63);
    CHECK(awc.walkers_per_crowd.size() == 4);
    CHECK(awc.walkers_per_rank[0] == 16);
    CHECK(awc.walkers_per_rank[3] == 15);
    CHECK(awc.walkers_per_crowd[0] == 4);
    CHECK(awc.walkers_per_crowd[3] == 3);

     awc = adjustGlobalWalkerCount(4,3,0,32,1.0,4);
    CHECK(awc.global_walkers == 128);
    CHECK(awc.walkers_per_crowd.size() == 4);
    CHECK(awc.walkers_per_rank[0] == 32);
    CHECK(awc.walkers_per_rank[3] == 32);
    CHECK(awc.walkers_per_crowd[0] == 8);
    CHECK(awc.walkers_per_crowd[3] == 8);
  }
  
  bool run() { return false; }

  int get_num_crowds() { return crowds_.size(); }
  
private:
  
};

}
} // namespace qmcplusplus
#endif
