//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_VMC_RENY_PARTICLEBYPARTICLE_UPDATE_H
#define QMCPLUSPLUS_VMC_RENY_PARTICLEBYPARTICLE_UPDATE_H
#include "QMCDrivers/QMCUpdateBase.h"

namespace qmcplusplus
{

/** @ingroup QMCDrivers  ParticleByParticle
*@brief Implements the VMC algorithm using particle-by-particle move with the drift equation.
*/
class VMCUpdateRenyiWithDriftFast: public QMCUpdateBase
{
public:
  /// Constructor.
  VMCUpdateRenyiWithDriftFast(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                              QMCHamiltonian& h, RandomGenerator_t& rg);

  ~VMCUpdateRenyiWithDriftFast();

  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);

private:
  std::vector<NewTimer*> myTimers;
};

}
#endif
