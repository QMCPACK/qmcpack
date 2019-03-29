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
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_DMC_UPDATE_PARTICLEBYPARTCLE_VMC_H
#define QMCPLUSPLUS_DMC_UPDATE_PARTICLEBYPARTCLE_VMC_H
#include "QMCDrivers/QMCUpdateBase.h"
#include "Utilities/NewTimer.h"

namespace qmcplusplus
{
class DMCUpdatePbyPVMC : public QMCUpdateBase
{
public:
  /// Constructor.
  DMCUpdatePbyPVMC(MCWalkerConfiguration& w,
                                 TrialWaveFunction& psi,
                                 QMCHamiltonian& h,
                                 RandomGenerator_t& rg);
  ///destructor
  ~DMCUpdatePbyPVMC();

  void advanceWalker(Walker_t& thisWalker, bool recompute);

private:
  TimerList_t myTimers;
};


} // namespace qmcplusplus

#endif
