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
#include "ContextForSteps.h"

namespace qmcplusplus
{

template<bool spinor>
ContextForSteps<spinor>::ContextForSteps(int num_walkers,
                                 int num_particles,
                                 std::vector<std::pair<int, int>> particle_group_indexes,
                                 RandomGenerator& random_gen)
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

  walker_deltas_.rs.resize(num_walkers * num_particles);
}

template class ContextForSteps<true>;
template class ContextForSteps<false>;

} // namespace qmcplusplus
