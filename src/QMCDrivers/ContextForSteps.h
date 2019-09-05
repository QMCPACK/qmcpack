//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 developers.
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

class MCPopulation;
struct DistanceTableData;
/** Thread local context for moving walkers
 *
 *  created once per driver per crowd
 *  might be merged with Crowd
 */
class ContextForSteps
{
public:
  using ParticlePositions = PtclOnLatticeTraits::ParticlePos_t;
  using MCPWalker  = Walker<QMCTraits, PtclOnLatticeTraits>;
  using RealType = QMCTraits::RealType;

  ContextForSteps(int num_walkers,
              int num_particles,
              std::vector<std::pair<int,int>> particle_group_indexes,
              RandomGenerator_t& random_gen);
  
  void loadCrowd(Crowd& crowd);

  std::vector<std::unique_ptr<ParticlePositions>>& get_walker_positions() { return walker_positions_; }

  int get_num_groups() const { return particle_group_indexes_.size(); }

  void nextDeltaRs() {
      std::for_each(walker_deltas_.begin(),
                walker_deltas_.end(),
                [this](auto& deltaR) {
                        makeGaussRandomWithEngine(*deltaR, random_gen_);
                });
  }
  
  std::vector<std::unique_ptr<ParticlePositions>>& get_walker_deltas() { return walker_deltas_; }

  int getPtclGroupStart(int group) const { return particle_group_indexes_[group].first; }
  int getPtclGroupEnd(int group) const { return particle_group_indexes_[group].first; }
  
protected:
// only one of these should exist
  ///Positions
  std::vector<std::unique_ptr<ParticlePositions>> walker_positions_;

  ///SoA copy of R
  std::vector<std::unique_ptr<VectorSoaContainer<RealType, OHMMS_DIM>>> positions_soa_;

  std::vector<std::unique_ptr<ParticlePositions>> walker_deltas_;
   
  /** indexes of start and stop of each particle group;
   *
   *  Seems like these should be iterators but haven't thought through the implications.
   */
  std::vector<std::pair<int,int>> particle_group_indexes_;

  /** map to handle distance tables
   *
   * myDistTableMap[source-particle-tag]= locator in the distance table
   * myDistTableMap[ObjectTag] === 0
   */
  std::map<std::string, int> myDistTableMap;

  /// distance tables that need to be updated by moving this ParticleSet
  std::vector<DistanceTableData*> DistTables;

  /// Descriptions from distance table creation.  Same order as DistTables.
  std::vector<std::string> distTableDescriptions;

  RandomGenerator_t& random_gen_;
};

}
#endif
