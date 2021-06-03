//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SOVMC_UPDATEALL_H
#define QMCPLUSPLUS_SOVMC_UPDATEALL_H
#include "QMCDrivers/QMCUpdateBase.h"

namespace qmcplusplus
{
/** @ingroup QMCDrivers  ParticleByParticle
 *@brief Implements the Spin VMC algorithm using particle-by-particle move.
 */
class SOVMCUpdateAll : public QMCUpdateBase
{
public:
  /// Constructor.
  SOVMCUpdateAll(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg);

  ~SOVMCUpdateAll();

  void advanceWalker(Walker_t& thisWalker, bool recompute);

private:
  /// Copy Constructor (disabled)
  SOVMCUpdateAll(const SOVMCUpdateAll&) = delete;
  /// Copy operator (disabled).
  SOVMCUpdateAll& operator=(const SOVMCUpdateAll&) = delete;
};

} // namespace qmcplusplus

#endif
