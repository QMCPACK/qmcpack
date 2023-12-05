//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020  QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

/** @file MCSample.h
 * @brief Stores particle configurations for later use in DMC and wavefunction optimization
 *
 * Stores less temporary data than the buffer object.
 */

#ifndef QMCPLUSPLUS_MCSAMPLE_H
#define QMCPLUSPLUS_MCSAMPLE_H

#include "Particle/ParticleSet.h"
#include "Particle/Walker.h"

namespace qmcplusplus
{
/** store minimum Walker data
 */
struct MCSample
{
  using WP       = WalkerProperties::Indexes;
  using Walker_t = ParticleSet::Walker_t;

  ParticleSet::ParticlePos R;
  ParticleSet::FullPrecRealType Weight;
  ParticleSet::ParticleScalar spins;
  ParticleSet::ParticleGradient G;
  ParticleSet::ParticleLaplacian L;
  ParticleSet::RealType LogPsi, Sign, PE, KE;

  inline MCSample(const ParticleSet& pset) : R(pset.R), spins(pset.spins) {}

  /// deprecated. Beyond w.R and w.spins, others are used perhaps somewhere but intended not to.
  inline MCSample(const Walker_t& w) : R(w.R), Weight(w.Weight), spins(w.spins), G(w.G), L(w.L)
  {
    LogPsi = w.Properties(WP::LOGPSI);
    Sign   = w.Properties(WP::SIGN);
    PE     = w.Properties(WP::LOCALPOTENTIAL);
    KE     = w.Properties(WP::LOCALENERGY) - PE;
  }

  inline size_t getNumPtcls() const { return R.size(); }

  inline MCSample(int n)
  {
    R.resize(n);
    spins.resize(n);
    G.resize(n);
    L.resize(n);
  }

  inline void convertToWalker(Walker_t& w) const
  {
    w.R                              = R;
    w.Weight                         = Weight;
    w.spins                          = spins;
    w.G                              = G;
    w.L                              = L;
    w.Properties(WP::LOGPSI)         = LogPsi;
    w.Properties(WP::SIGN)           = Sign;
    w.Properties(WP::LOCALPOTENTIAL) = PE;
    w.Properties(WP::LOCALENERGY)    = PE + KE;
  }
};

} // namespace qmcplusplus
#endif
