//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "MinimalParticlePool.h"
#include "Particle/VirtualParticleSet.h"


#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
TEST_CASE("VirtualParticleSet", "[particle]")
{
  auto pset_pool = MinimalParticlePool::make_NiO_a4(OHMMS::Controller);

  auto& ions  = *pset_pool.getParticleSet("i");
  auto& elecs = *pset_pool.getParticleSet("e");

  ions.addTable(ions);
  ions.update();
  elecs.addTable(ions);
  elecs.addTable(elecs);
  elecs.update();
}
} // namespace qmcplusplus
