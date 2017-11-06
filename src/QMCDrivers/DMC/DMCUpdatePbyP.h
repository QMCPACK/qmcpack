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
    
    


#ifndef QMCPLUSPLUS_DMC_UPDATE_PARTICLEBYPARTCLE_H
#define QMCPLUSPLUS_DMC_UPDATE_PARTICLEBYPARTCLE_H
#include "QMCDrivers/QMCUpdateBase.h"
#include "Utilities/NewTimer.h"
namespace qmcplusplus
{

class DMCUpdatePbyPWithRejectionFast: public QMCUpdateBase
{

public:

  /// Constructor.
  DMCUpdatePbyPWithRejectionFast(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                                 QMCHamiltonian& h, RandomGenerator_t& rg);
  ///destructor
  ~DMCUpdatePbyPWithRejectionFast();

  void advanceWalker(Walker_t& thisWalker, bool recompute);

private:
  TimerList_t myTimers;
};


enum DMCTimers
{
  DMC_buffer,
  DMC_movePbyP,
  DMC_hamiltonian,
  DMC_collectables,
  DMC_tmoves
};

extern TimerNameList_t<DMCTimers> DMCTimerNames;


}

#endif
