//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from ParticleSet.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include <type_traits>
#include "ContextForSteps.h"

namespace qmcplusplus
{
ContextForSteps::ContextForSteps(int num_walkers,
                                 int num_particles,
                                 std::vector<std::pair<int, int>> particle_group_indexes,
                                 RandomGenerator& random_gen)
    : particle_group_indexes_(particle_group_indexes), random_gen_(random_gen)
{}

} // namespace qmcplusplus
