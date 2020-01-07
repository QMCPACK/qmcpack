//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SODMC_UPDATE_PARTICLEBYPARTCLE_H
#define QMCPLUSPLUS_SODMC_UPDATE_PARTICLEBYPARTCLE_H
#include "QMCDrivers/QMCUpdateBase.h"
#include "Utilities/NewTimer.h"
namespace qmcplusplus
{
class SODMCUpdatePbyPWithRejectionFast : public QMCUpdateBase
{
public:
  /// Constructor.
  SODMCUpdatePbyPWithRejectionFast(MCWalkerConfiguration& w,
                                   TrialWaveFunction& psi,
                                   QMCHamiltonian& h,
                                   RandomGenerator_t& rg);
  ///destructor
  ~SODMCUpdatePbyPWithRejectionFast();

  void advanceWalker(Walker_t& thisWalker, bool recompute);

private:
  TimerList_t myTimers;
};


enum SODMCTimers
{
  SODMC_buffer,
  SODMC_movePbyP,
  SODMC_hamiltonian,
  SODMC_collectables,
  SODMC_tmoves
};

extern TimerNameList_t<SODMCTimers> SODMCTimerNames;


} // namespace qmcplusplus

#endif
