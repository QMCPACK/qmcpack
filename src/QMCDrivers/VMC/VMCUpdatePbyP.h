//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_VMC_PARTICLEBYPARTICLE_UPDATE_H
#define QMCPLUSPLUS_VMC_PARTICLEBYPARTICLE_UPDATE_H
#include "QMCDrivers/QMCUpdateBase.h"

namespace qmcplusplus
{

void add_vmc_timers(std::vector<NewTimer*>& timers);

/** @ingroup QMCDrivers  ParticleByParticle
 *@brief Implements the VMC algorithm using particle-by-particle move.
 */
class VMCUpdatePbyP: public QMCUpdateBase
{
public:
  /// Constructor.
  VMCUpdatePbyP(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                QMCHamiltonian& h, RandomGenerator_t& rg);

  ~VMCUpdatePbyP();

  void advanceWalker(Walker_t& thisWalker, bool recompute);

  RealType advanceWalkerForEE(Walker_t& w1, std::vector<PosType>& dR, std::vector<int>& iats, std::vector<int>& rs, std::vector<RealType>& ratios);

private:
  std::vector<NewTimer*> myTimers;
};

/** @ingroup QMCDrivers  ParticleByParticle
 *@brief Implements the VMC algorithm using particle-by-particle move with the drift equation.
 */
class VMCUpdatePbyPWithDriftFast: public QMCUpdateBase
{
public:
  /// Constructor.
  VMCUpdatePbyPWithDriftFast(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                             QMCHamiltonian& h, RandomGenerator_t& rg);

  ~VMCUpdatePbyPWithDriftFast();

  void advanceWalker(Walker_t& thisWalker, bool recompute);

  RealType advanceWalkerForEE(Walker_t& w1, std::vector<PosType>& dR, std::vector<int>& iats, std::vector<int>& rs, std::vector<RealType>& ratios);
private:
  std::vector<NewTimer*> myTimers;
};

}

#endif
