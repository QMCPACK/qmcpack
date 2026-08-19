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

#include <stdexcept>

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
void PSdispatcher::flex_makeMoveAllParticles(const RefVectorWithLeader<ParticleSet>& p_list,
                                             const RefVector<Walker_t>& resolved_walkers,
                                             const MCCoords<CT>& displacements,
                                             std::vector<bool>& are_valid,
                                             bool skipSK) const
{
  if (use_batch_)
  {
    ParticleSet::mw_makeMoveAllParticles(p_list, resolved_walkers, displacements, are_valid, skipSK);
    return;
  }

  const size_t num_walkers = p_list.size();
  if (num_walkers == 0 || resolved_walkers.size() != num_walkers || are_valid.size() != num_walkers)
    throw std::runtime_error("All-particle transaction walker counts do not match.");
  const size_t num_particles = p_list[0].getTotalNum();
  if (displacements.positions.size() != num_particles * num_walkers)
    throw std::runtime_error("Flattened all-particle displacement count does not match the crowd.");
  if constexpr (CT == CoordsType::POS_SPIN)
    if (displacements.spins.size() != num_particles * num_walkers)
      throw std::runtime_error("Flattened all-particle spin displacement count does not match the crowd.");

  // Complete crowd validation prevents a malformed later walker from exposing a partial update.
  for (size_t iw = 0; iw < num_walkers; ++iw)
  {
    const Walker_t& resolved_walker = resolved_walkers[iw];
    if (p_list[iw].getActivePtcl() != -1)
      throw std::runtime_error("Cannot start an all-particle transaction with an active particle.");
    if (p_list[iw].getTotalNum() != num_particles || p_list[iw].R.size() != num_particles ||
        p_list[iw].spins.size() != num_particles || resolved_walker.R.size() != num_particles ||
        resolved_walker.spins.size() != num_particles)
      throw std::runtime_error("All-particle transaction particle counts do not match.");
  }

  for (size_t iw = 0; iw < num_walkers; ++iw)
  {
    MCCoords<CT> walker_displacements(num_particles);
    for (size_t iat = 0; iat < num_particles; ++iat)
    {
      const size_t flat_index             = iat * num_walkers + iw;
      walker_displacements.positions[iat] = displacements.positions[flat_index];
      if constexpr (CT == CoordsType::POS_SPIN)
        walker_displacements.spins[iat] = displacements.spins[flat_index];
    }
    are_valid[iw] = p_list[iw].makeMoveAllParticles(resolved_walkers[iw], walker_displacements, skipSK);
  }
}

void PSdispatcher::flex_accept_rejectMoveAllParticles(const RefVectorWithLeader<ParticleSet>& p_list,
                                                      const RefVector<Walker_t>& resolved_walkers,
                                                      const std::vector<bool>& accepted,
                                                      bool skipSK) const
{
  if (use_batch_)
  {
    ParticleSet::mw_accept_rejectMoveAllParticles(p_list, resolved_walkers, accepted, skipSK);
    return;
  }

  const size_t num_walkers = p_list.size();
  if (num_walkers == 0 || resolved_walkers.size() != num_walkers || accepted.size() != num_walkers)
    throw std::runtime_error("All-particle transaction walker counts do not match.");
  const size_t num_particles = p_list[0].getTotalNum();
  for (size_t iw = 0; iw < num_walkers; ++iw)
  {
    const Walker_t& resolved_walker = resolved_walkers[iw];
    if (p_list[iw].getActivePtcl() != -1)
      throw std::runtime_error("Cannot resolve an all-particle transaction with an active particle.");
    if (p_list[iw].getTotalNum() != num_particles || p_list[iw].R.size() != num_particles ||
        p_list[iw].spins.size() != num_particles || resolved_walker.R.size() != num_particles ||
        resolved_walker.spins.size() != num_particles)
      throw std::runtime_error("All-particle transaction particle counts do not match.");
  }

  for (size_t iw = 0; iw < num_walkers; ++iw)
    p_list[iw].accept_rejectMoveAllParticles(resolved_walkers[iw], accepted[iw], skipSK);
}

template<CoordsType CT>
void PSdispatcher::flex_makeMove(const RefVectorWithLeader<ParticleSet>& p_list,
                                 int iat,
                                 const MCCoords<CT>& displs,
                                 std::vector<bool>& are_valid) const
{
  if (use_batch_)
    ParticleSet::mw_makeMove(p_list, iat, displs, are_valid);
  else
    for (size_t iw = 0; iw < p_list.size(); iw++)
      if constexpr (CT == CoordsType::POS_SPIN)
        are_valid[iw] = p_list[iw].makeMoveAndCheckWithSpin(iat, displs.positions[iw], displs.spins[iw]);
      else
        are_valid[iw] = p_list[iw].makeMoveAndCheck(iat, displs.positions[iw]);
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

template void PSdispatcher::flex_makeMoveAllParticles<CoordsType::POS>(const RefVectorWithLeader<ParticleSet>& p_list,
                                                                       const RefVector<Walker_t>& resolved_walkers,
                                                                       const MCCoords<CoordsType::POS>& displacements,
                                                                       std::vector<bool>& are_valid,
                                                                       bool skipSK) const;
template void PSdispatcher::flex_makeMoveAllParticles<CoordsType::POS_SPIN>(
    const RefVectorWithLeader<ParticleSet>& p_list,
    const RefVector<Walker_t>& resolved_walkers,
    const MCCoords<CoordsType::POS_SPIN>& displacements,
    std::vector<bool>& are_valid,
    bool skipSK) const;
template void PSdispatcher::flex_makeMove<CoordsType::POS>(const RefVectorWithLeader<ParticleSet>& p_list,
                                                           int iat,
                                                           const MCCoords<CoordsType::POS>& displs,
                                                           std::vector<bool>& are_valid) const;
template void PSdispatcher::flex_makeMove<CoordsType::POS_SPIN>(const RefVectorWithLeader<ParticleSet>& p_list,
                                                                int iat,
                                                                const MCCoords<CoordsType::POS_SPIN>& displs,
                                                                std::vector<bool>& are_valid) const;
template void PSdispatcher::flex_accept_rejectMove<CoordsType::POS>(const RefVectorWithLeader<ParticleSet>& p_list,
                                                                    int iat,
                                                                    const std::vector<bool>& isAccepted,
                                                                    bool forward_mode) const;
template void PSdispatcher::flex_accept_rejectMove<CoordsType::POS_SPIN>(const RefVectorWithLeader<ParticleSet>& p_list,
                                                                         int iat,
                                                                         const std::vector<bool>& isAccepted,
                                                                         bool forward_mode) const;
} // namespace qmcplusplus
