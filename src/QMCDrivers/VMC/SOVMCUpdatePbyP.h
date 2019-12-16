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
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SOVMC_PARTICLEBYPARTICLE_UPDATE_H
#define QMCPLUSPLUS_SOVMC_PARTICLEBYPARTICLE_UPDATE_H
#include "QMCDrivers/QMCUpdateBase.h"

namespace qmcplusplus
{
/** @ingroup QMCDrivers  ParticleByParticle
 *@brief Implements the VMC algorithm using particle-by-particle move, including spin-moves
 */
class SOVMCUpdatePbyP : public QMCUpdateBase
{
public:
  /// Constructor.
  SOVMCUpdatePbyP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg);

  ~SOVMCUpdatePbyP();

  void advanceWalker(Walker_t& thisWalker, bool recompute);

private:
  NewTimer& buffer_timer_;
  NewTimer& movepbyp_timer_;
  NewTimer& hamiltonian_timer_;
  NewTimer& collectables_timer_;
};

} // namespace qmcplusplus

#endif
