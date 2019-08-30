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
#include "QMCDrivers/ContextForSteps.h"
#include "QMCDrivers/MCPopulation.h"

namespace qmcplusplus
{
ContextForSteps::ContextForSteps(int num_walkers,
                         int num_particles,
                         std::vector<std::pair<int, int>> particle_group_indexes,
                         RandomGenerator_t& random_gen)
    : particle_group_indexes_(particle_group_indexes), random_gen_(random_gen)
{
  /** glambda to create type T with constructor T(int) and put in it unique_ptr
   *
   *  captures num_particles to use as argument to constructor
   *  gets T for type unique_ptr unique is templated on
   */
  auto constructT = [num_particles](auto& unique) {
    unique.reset(new typename std::remove_pointer<decltype(unique.get())>::type(num_particles));
  };

  walker_positions_.resize(num_walkers);
  std::for_each(walker_positions_.begin(), walker_positions_.end(), constructT);
  walker_deltas_.resize(num_walkers);
  std::for_each(walker_deltas_.begin(), walker_deltas_.end(), constructT);
  positions_soa_.resize(num_walkers);
  std::for_each(positions_soa_.begin(), positions_soa_.end(), constructT);
}

void ContextForSteps::loadCrowd(Crowd& crowd)
{
  auto it_walker        = crowd.beginWalkers();
  auto it_positions     = walker_positions_.begin();
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
