//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from ParticleSet.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include <type_traits>
#include "Particle/MoveContext.h"
#include "Particle/MCPopulation.h"

namespace qmcplusplus
{
MoveContext::MoveContext(Crowd& crowd, int particles, std::vector<std::pair<int, int>> particle_group_indexes, RandomGenerator_t& random_gen)
    : particle_group_indexes_(particle_group_indexes), random_gen_(random_gen)
{
  // std::for_each(positions_per_walker_.begin(), positions_per_walker_.end(),
  //               [particles](std::unique_ptr<ParticlePositions>& positions_) {
  //                 positions_ = std::make_unique<ParticlePositions>(particles);
  //               });

  auto vecUniqueResize = [particles](auto& vec_unique) {
                           vec_unique.reset(new typename std::remove_pointer<decltype(vec_unique.get())>::type (particles));
  };

  positions_per_walker_.resize(crowd.size());
  std::for_each(positions_per_walker_.begin(), positions_per_walker_.end(), vecUniqueResize);
  delta_positions_per_walker_.resize(crowd.size());
  std::for_each(delta_positions_per_walker_.begin(), delta_positions_per_walker_.end(), vecUniqueResize);
  
  positions_soa_.resize(crowd.size());
  std::for_each(positions_soa_.begin(), positions_soa_.end(), vecUniqueResize);
  
  // std::for_each(positions_soa_.begin(), positions_soa_.end(),
  //               [particles](std::unique_ptr<VectorSoaContainer<RealType, OHMMS_DIM>>& walker_positions_soa) {
  //                 walker_positions_soa = std::make_unique<VectorSoaContainer<RealType, OHMMS_DIM>>(particles);
  //               });

  
  
}

void MoveContext::loadCrowd(Crowd& crowd)
{
  auto it_walker        = crowd.beginWalkers();
  auto it_positions     = positions_per_walker_.begin();
  auto it_positions_soa = positions_soa_.begin();
  while (it_walker != crowd.endWalkers())
  {
    **it_positions = it_walker->get().R;
    (*it_positions_soa)->copyIn(**it_positions);
    ++it_walker;
    ++it_positions;
    ++it_positions_soa;
  }
  //positions_soa_.copyIn(positions_);
  // in certain cases, full tables must be ready
  // for (int i = 0; i < DistTables.size(); i++)
  //     if (DistTables[i]->DTType == DT_AOS || DistTables[i]->Need_full_table_loadWalker)
  //       DistTables[i]->evaluate(*this);
  //   //computed so that other objects can use them, e.g., kSpaceJastrow
  //   if (SK && SK->DoUpdate)
  //     SK->UpdateAllPart(*this);

  // activePtcl = -1;
}

} // namespace qmcplusplus
