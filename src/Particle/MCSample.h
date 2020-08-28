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

  ParticleSet::ParticlePos_t R;
  ParticleSet::ParticleScalar_t spins;
  ParticleSet::ParticleGradient_t G;
  ParticleSet::ParticleLaplacian_t L;
  ParticleSet::RealType LogPsi, Sign, PE, KE;

  inline MCSample(const Walker_t& w) : R(w.R), spins(w.spins), G(w.G), L(w.L)
  {
    LogPsi = w.Properties(WP::LOGPSI);
    Sign   = w.Properties(WP::SIGN);
    PE     = w.Properties(WP::LOCALPOTENTIAL);
    KE     = w.Properties(WP::LOCALENERGY) - PE;
  }

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
