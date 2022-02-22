//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from ParticleSet.h
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_MOVECONTEXT_H
#define QMCPLUSPLUS_MOVECONTEXT_H

#include <map>
#include <vector>
#include "OhmmsSoA/VectorSoaContainer.h"
#include "Configuration.h"
#include "Particle/Walker.h"
#include "QMCDrivers/Crowd.h"
#include "ParticleBase/RandomSeqGenerator.h"

namespace qmcplusplus
{
/** Thread local context for moving walkers
 *
 *  created once per driver per crowd
 *  It's two significant responsibilities are holding the thread local RandomGen_t
 *  And the particle group indexes.
 *
 *  
 */
class ContextForSteps
{
public:
  using ParticlePositions = PtclOnLatticeTraits::ParticlePos;
  using PosType           = QMCTraits::PosType;
  using MCPWalker         = Walker<QMCTraits, PtclOnLatticeTraits>;
  using RealType          = QMCTraits::RealType;

  ContextForSteps(int num_walkers,
                  int num_particles,
                  std::vector<std::pair<int, int>> particle_group_indexes,
                  RandomGenerator& random_gen);

  int get_num_groups() const { return particle_group_indexes_.size(); }
  RandomGenerator& get_random_gen() { return random_gen_; }

  int getPtclGroupStart(int group) const { return particle_group_indexes_[group].first; }
  int getPtclGroupEnd(int group) const { return particle_group_indexes_[group].second; }

protected:
  /** indexes of start and stop of each particle group;
   *
   *  Seems like these should be iterators but haven't thought through the implications.
   */
  std::vector<std::pair<int, int>> particle_group_indexes_;

  RandomGenerator& random_gen_;
};

} // namespace qmcplusplus
#endif
