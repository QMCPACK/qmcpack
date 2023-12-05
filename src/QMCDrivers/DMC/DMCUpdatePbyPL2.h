//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov,  Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_DMC_UPDATE_PARTICLEBYPARTCLE_L2_H
#define QMCPLUSPLUS_DMC_UPDATE_PARTICLEBYPARTCLE_L2_H
#include "QMCDrivers/QMCUpdateBase.h"
#include "Utilities/TimerManager.h"

namespace qmcplusplus
{
class DMCUpdatePbyPL2 : public QMCUpdateBase
{
public:
  /// Constructor.
  DMCUpdatePbyPL2(MCWalkerConfiguration& w,
                  TrialWaveFunction& psi,
                  QMCHamiltonian& h,
                  RandomBase<FullPrecRealType>& rg);
  ///destructor
  ~DMCUpdatePbyPL2() override;

  void advanceWalker(Walker_t& thisWalker, bool recompute) override;

private:
  TimerList_t myTimers;
};

//extern TimerNameList_t<DMCTimers> DMCTimerNames;


} // namespace qmcplusplus

#endif
