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
  SOVMCUpdatePbyP(MCWalkerConfiguration& w,
                  TrialWaveFunction& psi,
                  QMCHamiltonian& h,
                  RandomBase<FullPrecRealType>& rg);

  ~SOVMCUpdatePbyP() override;

  void advanceWalker(Walker_t& thisWalker, bool recompute) override;

private:
  NewTimer& buffer_timer_;
  NewTimer& movepbyp_timer_;
  NewTimer& hamiltonian_timer_;
  NewTimer& collectables_timer_;
};

} // namespace qmcplusplus

#endif
