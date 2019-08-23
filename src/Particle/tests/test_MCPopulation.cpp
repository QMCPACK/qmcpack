//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"
#include "Configuration.h"
#include "OhmmsPETE/TinyVector.h"
#include "Particle/MCPopulation.h"

namespace qmcplusplus
{
TEST_CASE("MCPopulation::createWalkers", "[particle][population]")
{
  MCPopulation population(4, 2);
  ParticleAttrib<TinyVector<QMCTraits::RealType, 3>> some_pos(2);
  some_pos[0] = TinyVector<double, 3>(1.0, 0.0, 0.0);
  some_pos[1] = TinyVector<double, 3>(0.0, 1.0, 0.0);

  population.createWalkers(8, some_pos);
  REQUIRE(population.get_walkers().size() == 8);
}

TEST_CASE("MCPopulation::createWalkers first touch", "[particle][population]")
{
  MCPopulation population(1, 2);
  ParticleAttrib<TinyVector<QMCTraits::RealType, 3>> some_pos(2);
  some_pos[0] = TinyVector<double, 3>(1.0, 0.0, 0.0);
  some_pos[1] = TinyVector<double, 3>(0.0, 1.0, 0.0);

  population.createWalkers(8, 2, 16, some_pos);
  REQUIRE(population.get_walkers().size() == 16);
  // Someday here we should use hwloc or similar to query that the memory is close to
  // each thread.
}


} // namespace qmcplusplus
