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

  /** Pass an all particle move directly to the multiwalker implementation.
   *
   * @tparam CT coordinate type, either POS or POS_SPIN
   * @param[in,out] p_list ParticleSets to move
   * @param[in] displacements particle displacements in particle-major order
   * @param[out] are_valid lattice-validity result for each particle move, in particle-major order
   * @throws std::runtime_error if this dispatcher is not configured for batched execution
   */
  template<CoordsType CT>
  void flex_makeMoveAllParticles(const RefVectorWithLeader<ParticleSet>& p_list,
                                 const MCCoords<CT>& displacements,
                                 std::vector<bool>& are_valid) const;

  /** Pass all particle move decisions directly to the multiwalker implementation.
   *
   * @param[in,out] p_list ParticleSets containing the attempted moves
   * @param[in] accepted one acceptance decision for each ParticleSet
   * @throws std::runtime_error if this dispatcher is not configured for batched execution
   */
  void flex_accept_rejectMoveAllParticles(const RefVectorWithLeader<ParticleSet>& p_list,
                                          const std::vector<bool>& accepted) const;

  template<CoordsType CT>
  void flex_makeMove(const RefVectorWithLeader<ParticleSet>& p_list,
                     int iat,
                     const MCCoords<CT>& displs,
                     std::vector<bool>& are_valid) const;

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
