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
#include "ResourceCollection.h"
#include "DistanceTable.h"

namespace qmcplusplus
{
TEST_CASE("VirtualParticleSet", "[particle]")
{
  auto pset_pool = MinimalParticlePool::make_NiO_a4(OHMMS::Controller);

  auto& ions  = *pset_pool.getParticleSet("i");
  auto& elecs = *pset_pool.getParticleSet("e");

  elecs.R[0] = {1, 2, 3};
  elecs.R[1] = {2, 1, 3};
  elecs.R[2] = {3, 1, 2};
  elecs.R[3] = {3, 2, 1};

  ions.addTable(ions);
  ions.update();
  elecs.addTable(ions);
  elecs.addTable(elecs);
  elecs.update();


  ParticleSet elecs_clone(elecs);
  elecs_clone.update();

  VirtualParticleSet vp_Ni(elecs, 3);
  VirtualParticleSet vp_O(elecs, 5);

  VirtualParticleSet vp_Ni_clone(elecs_clone, 3);
  VirtualParticleSet vp_O_clone(elecs_clone, 5);

  vp_Ni_clone.makeMoves(elecs_clone, 3, {{0.1, 0.2, 0.3}, {0.2, 0.1, 0.3}, {0.3, 0.1, 0.2}});
  const DistanceTableAB& dt_vp_ion = vp_Ni_clone.getDistTableAB(0);
  CHECK(Approx(dt_vp_ion.getDistances()[2][1]) == 2.5020600118);

  // two walkers form a workgroup.
  // One electron of the first walker gets inside O sphere and one electron of the other gets inside Ni sphere.
  RefVectorWithLeader<VirtualParticleSet> vp_list(vp_Ni, {vp_O, vp_Ni_clone});
  ResourceCollection collection{"NLPPcollection"};
  vp_Ni.createResource(collection);

  {
    ResourceCollectionTeamLock<VirtualParticleSet> vp_res_lock(collection, vp_list);
  }

  vp_Ni_clone.makeMoves(elecs_clone, 3, {{0.1, 0.2, 0.3}, {0.3, 0.1, 0.2}, {0.2, 0.1, 0.3}});
  CHECK(Approx(vp_Ni_clone.R[2][0]) == 3.2);
  CHECK(Approx(vp_Ni_clone.R[2][1]) == 2.1);
  CHECK(Approx(vp_Ni_clone.R[2][2]) == 1.3);

  REQUIRE(dt_vp_ion.getDistances().size() == 3);
  CHECK(Approx(dt_vp_ion.getDistances()[2][1]) == 2.5784519198);
}
} // namespace qmcplusplus
