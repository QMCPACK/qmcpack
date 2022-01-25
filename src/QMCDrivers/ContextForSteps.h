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

/** Thread local context for moving walkers
 *
 *  created once per driver per crowd
 *  It's two significant responsibilities are holding the thread local RandomGen_t
 *  And the particle group indexes.
 *
 *  
 */
template<bool spinor = false>
class ContextForSteps
{
public:
  using ParticlePositions = PtclOnLatticeTraits::ParticlePos_t;
  using PosType           = QMCTraits::PosType;
  using MCPWalker         = Walker<QMCTraits, PtclOnLatticeTraits>;
  using RealType          = QMCTraits::RealType;

  enum class MCCoordsTypes
  {
    RS,
    RSSPINS
  };

  static constexpr MCCoordsTypes translated_ct = spinor ? MCCoordsTypes::RSSPINS : MCCoordsTypes::RS;

  template<MCCoordsTypes CT = MCCoordsTypes::RS>
  struct MCCoords
  {
    std::vector<QMCTraits::PosType> rs;
  };

  template<>
  struct MCCoords<MCCoordsTypes::RSSPINS>
  {
    std::vector<QMCTraits::PosType> rs;
    std::vector<std::complex<double>> spins;
  };

  template<MCCoordsTypes CT = MCCoordsTypes::RS>
  struct MCCIt
  {
    std::vector<QMCTraits::PosType>::iterator irs;
  };

  template<>
  struct MCCIt<MCCoordsTypes::RSSPINS>
  {
    std::vector<QMCTraits::PosType>::iterator irs;
    std::vector<std::complex<double>>::iterator spins;
  };

  ContextForSteps(int num_walkers,
                  int num_particles,
                  std::vector<std::pair<int, int>> particle_group_indexes,
                  RandomGenerator& random_gen);

  int get_num_groups() const { return particle_group_indexes_.size(); }
  RandomGenerator& get_random_gen() { return random_gen_; }

  void nextDeltas(size_t num_rs)
  {
    walker_deltas_.rs.resize(num_rs);
    makeGaussRandomWithEngine(walker_deltas_.rs, random_gen_);
    // hate to repeat this pattern, this should never resize.
    if constexpr (std::is_same<decltype(walker_deltas_), MCCoords<MCCoordsTypes::RSSPINS>>::value)
    {
      walker_deltas_.spins.resize(num_rs);
      makeGaussRandomWithEngine(walker_deltas_.spins, random_gen_);
    }
  }

  MCCoords<translated_ct>& get_walker_deltas() { return walker_deltas_; }

  MCCIt<translated_ct> deltasBegin()
  {
    if constexpr (std::is_same<decltype(walker_deltas_), MCCoords<MCCoordsTypes::RS>>::value)
                     return {walker_deltas_.rs.begin()};
    else
        return {walker_deltas_.rs.begin(), walker_deltas_.spins.begin()};

  };

  int getPtclGroupStart(int group) const { return particle_group_indexes_[group].first; }
  int getPtclGroupEnd(int group) const { return particle_group_indexes_[group].second; }

protected:
  MCCoords<translated_ct> walker_deltas_;

  /** indexes of start and stop of each particle group;
   *
   *  Seems like these should be iterators but haven't thought through the implications.
   */
  std::vector<std::pair<int, int>> particle_group_indexes_;

  RandomGenerator& random_gen_;
};

extern template class ContextForSteps<true>;
extern template class ContextForSteps<false>;

} // namespace qmcplusplus
#endif
