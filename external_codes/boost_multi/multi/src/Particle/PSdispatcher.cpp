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
                                   const std::vector<bool>& recompute,
                                   bool pbyp) const
{
  if (use_batch_)
    ParticleSet::mw_loadWalker(p_list, walkers, recompute, pbyp);
  else
    for (size_t iw = 0; iw < p_list.size(); iw++)
      if (recompute[iw])
      {
        p_list[iw].loadWalker(walkers[iw], false);
        // loadWalker and mw_loadWalker doesn't have the same behavior. Need the following update call.
        p_list[iw].update(true);
      }
}

void PSdispatcher::flex_update(const RefVectorWithLeader<ParticleSet>& p_list, bool skipSK) const
{
  if (use_batch_)
    ParticleSet::mw_update(p_list, skipSK);
  else
    for (ParticleSet& pset : p_list)
      pset.update(skipSK);
}

template<CoordsType CT>
void PSdispatcher::flex_makeMove(const RefVectorWithLeader<ParticleSet>& p_list,
                                 int iat,
                                 const MCCoords<CT>& displs) const
{
  if (use_batch_)
    ParticleSet::mw_makeMove(p_list, iat, displs);
  else
    for (size_t iw = 0; iw < p_list.size(); iw++)
      if constexpr (CT == CoordsType::POS_SPIN)
        p_list[iw].makeMoveWithSpin(iat, displs.positions[iw], displs.spins[iw]);
      else
        p_list[iw].makeMove(iat, displs.positions[iw]);
}

template<CoordsType CT>
void PSdispatcher::flex_accept_rejectMove(const RefVectorWithLeader<ParticleSet>& p_list,
                                          int iat,
                                          const std::vector<bool>& isAccepted,
                                          bool forward_mode) const
{
  if (use_batch_)
    ParticleSet::mw_accept_rejectMove<CT>(p_list, iat, isAccepted, forward_mode);
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

template void PSdispatcher::flex_makeMove<CoordsType::POS>(const RefVectorWithLeader<ParticleSet>& p_list,
                                                           int iat,
                                                           const MCCoords<CoordsType::POS>& displs) const;
template void PSdispatcher::flex_makeMove<CoordsType::POS_SPIN>(const RefVectorWithLeader<ParticleSet>& p_list,
                                                                int iat,
                                                                const MCCoords<CoordsType::POS_SPIN>& displs) const;
template void PSdispatcher::flex_accept_rejectMove<CoordsType::POS>(const RefVectorWithLeader<ParticleSet>& p_list,
                                                                    int iat,
                                                                    const std::vector<bool>& isAccepted,
                                                                    bool forward_mode) const;
template void PSdispatcher::flex_accept_rejectMove<CoordsType::POS_SPIN>(const RefVectorWithLeader<ParticleSet>& p_list,
                                                                         int iat,
                                                                         const std::vector<bool>& isAccepted,
                                                                         bool forward_mode) const;
} // namespace qmcplusplus
