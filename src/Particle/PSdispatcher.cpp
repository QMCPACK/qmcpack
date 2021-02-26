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


#include "PSdispatcher.h"

namespace qmcplusplus
{
PSdispatcher::PSdispatcher(bool use_batch) : use_batch_(use_batch) {}

void PSdispatcher::flex_loadWalker(const RefVectorWithLeader<ParticleSet>& p_list,
                                   const RefVector<Walker_t>& walkers,
                                   bool pbyp) const
{
  if (use_batch_)
    ParticleSet::mw_loadWalker(p_list, walkers, pbyp);
  else
    for (size_t iw = 0; iw < p_list.size(); iw++)
      p_list[iw].loadWalker(walkers[iw], pbyp);
}

void PSdispatcher::flex_update(const RefVectorWithLeader<ParticleSet>& p_list, bool skipSK) const
{
  if (use_batch_)
    ParticleSet::mw_update(p_list, skipSK);
  else
    for (ParticleSet& pset : p_list)
      pset.update(skipSK);
}

void PSdispatcher::flex_makeMove(const RefVectorWithLeader<ParticleSet>& p_list,
                                 int iat,
                                 const std::vector<SingleParticlePos_t>& displs) const
{
  if (use_batch_)
    ParticleSet::mw_makeMove(p_list, iat, displs);
  else
    for (size_t iw = 0; iw < p_list.size(); iw++)
      p_list[iw].makeMove(iat, displs[iw]);
}

void PSdispatcher::flex_accept_rejectMove(const RefVectorWithLeader<ParticleSet>& p_list,
                                          int iat,
                                          const std::vector<bool>& isAccepted,
                                          bool forward_mode) const
{
  if (use_batch_)
    ParticleSet::mw_accept_rejectMove(p_list, iat, isAccepted, forward_mode);
  else
    for (size_t iw = 0; iw < p_list.size(); iw++)
      p_list[iw].accept_rejectMove(iat, isAccepted[iw], forward_mode);
}

void PSdispatcher::flex_donePbyP(const RefVectorWithLeader<ParticleSet>& p_list) const
{
  if (use_batch_)
    ParticleSet::mw_donePbyP(p_list);
  else
    for (ParticleSet& pset : p_list)
      pset.donePbyP();
}

void PSdispatcher::flex_saveWalker(const RefVectorWithLeader<ParticleSet>& p_list,
                                   const RefVector<Walker_t>& walkers) const
{
  if (use_batch_)
    ParticleSet::mw_saveWalker(p_list, walkers);
  else
    for (size_t iw = 0; iw < p_list.size(); iw++)
      p_list[iw].saveWalker(walkers[iw]);
}

} // namespace qmcplusplus
