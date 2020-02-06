//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, , doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_QMCDRIVERNEWTESTWRAPPER_H
#define QMCPLUSPLUS_QMCDRIVERNEWTESTWRAPPER_H
#include "QMCDrivers/QMCDriverNew.h"
#include "QMCDrivers/DriverTraits.h"

namespace qmcplusplus
{
/** Unit testing an impure virtual base class
 *  requires a absolute minimal subtype
 */
class QMCDriverNewTestWrapper : public QMCDriverNew
{
public:
  using Base = QMCDriverNew;
  QMCDriverNewTestWrapper(QMCDriverInput&& input,
                          MCPopulation& population,
                          TrialWaveFunction& psi,
                          QMCHamiltonian& h,
                          WaveFunctionPool& ppool,
                          Communicate* comm)
    : QMCDriverNew(std::move(input), population, psi, h, ppool, "QMCDriverTestWrapper::", comm)
  {}

  QMCRunType getRunType() { return QMCRunType::DUMMY; }

  void set_calc_walker_answers(int walkers_per_crowd, int local_walkers_fake)
    { walkers_per_crowd_fake_=walkers_per_crowd;
      local_walkers_fake_ = local_walkers_fake; }
  // Notice that this is a crap method in that we have to fake all the side effects of the
  // calculation in the child class.
  IndexType calc_default_local_walkers(IndexType walkers_per_rank)
  {
    walkers_per_crowd_ = walkers_per_crowd_fake_;
    return local_walkers_fake_;
  }

  void process(xmlNodePtr node)
  {
    IndexType local_walkers = calc_default_local_walkers(qmcdriver_input_.get_walkers_per_rank());
  
  // side effect updates walkers_per_crowd_;
    makeLocalWalkers(local_walkers, ParticleAttrib<TinyVector<QMCTraits::RealType, 3>>(population_.get_num_particles()));

    Base::process(node);
  }
  
  bool run() { return false; }
private:
  int walkers_per_crowd_fake_ = 4;
  int local_walkers_fake_ = 32;  
};

} // namespace qmcplusplus
#endif
