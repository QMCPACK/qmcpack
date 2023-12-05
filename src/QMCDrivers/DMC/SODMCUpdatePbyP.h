//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories 
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SODMC_UPDATE_PARTICLEBYPARTCLE_H
#define QMCPLUSPLUS_SODMC_UPDATE_PARTICLEBYPARTCLE_H
#include "QMCDrivers/QMCUpdateBase.h"
#include "Utilities/TimerManager.h"
namespace qmcplusplus
{
class SODMCUpdatePbyPWithRejectionFast : public QMCUpdateBase
{
public:
  /// Constructor.
  SODMCUpdatePbyPWithRejectionFast(MCWalkerConfiguration& w,
                                   TrialWaveFunction& psi,
                                   QMCHamiltonian& h,
                                   RandomBase<FullPrecRealType>& rg);
  ///destructor
  ~SODMCUpdatePbyPWithRejectionFast() override;

  void advanceWalker(Walker_t& thisWalker, bool recompute) override;

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

extern const TimerNameList_t<SODMCTimers> SODMCTimerNames;


} // namespace qmcplusplus

#endif
