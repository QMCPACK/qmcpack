//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_PSDISPATCH_H
#define QMCPLUSPLUS_PSDISPATCH_H

#include "ParticleSet.h"

namespace qmcplusplus
{
/** Wrappers for dispatching to ParticleSet single walker APIs or mw_ APIs.
 * This should be only used by QMC drivers.
 * member function names must match mw_ APIs in TrialWaveFunction
 */
class PSdispatcher
{
public:
  using Walker_t          = ParticleSet::Walker_t;
  using SingleParticlePos = ParticleSet::SingleParticlePos;
  using Scalar_t          = ParticleSet::Scalar_t;

  PSdispatcher(bool use_batch);

  void flex_loadWalker(const RefVectorWithLeader<ParticleSet>& p_list,
                       const RefVector<Walker_t>& walkers,
                       const std::vector<bool>& recompute,
                       bool pbyp) const;

  void flex_update(const RefVectorWithLeader<ParticleSet>& p_list, bool skipSK = false) const;

  template<CoordsType CT>
  void flex_makeMove(const RefVectorWithLeader<ParticleSet>& p_list, int iat, const MCCoords<CT>& displs) const;

  template<CoordsType CT>
  void flex_accept_rejectMove(const RefVectorWithLeader<ParticleSet>& p_list,
                              int iat,
                              const std::vector<bool>& isAccepted,
                              bool forward_mode = true) const;

  void flex_saveWalker(const RefVectorWithLeader<ParticleSet>& p_list, const RefVector<Walker_t>& walkers) const;

  void flex_donePbyP(const RefVectorWithLeader<ParticleSet>& p_list) const;

private:
  bool use_batch_;
};
} // namespace qmcplusplus

#endif
