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
#include "QMCDrivers/MoveContext.h"
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/TinyVectorOps.h"
#include "Estimators/tests/FakeEstimator.h"

namespace qmcplusplus
{
using WalkerMCP = Walker<QMCTraits, PtclOnLatticeTraits>;

TEST_CASE("MoveContext::loadWalker", "[particle]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate* comm = OHMMS::Controller;

  EstimatorManagerBase em(comm);
  FakeEstimator* fake_est = new FakeEstimator;
  em.add(fake_est, "fake");

  Crowd crowd(em);

  int num_particles = 1;
  WalkerMCP walker1(num_particles);
  walker1.R[0] = TinyVector<double, 3>(1.0, 0.0, 0.0);
  WalkerMCP walker2(num_particles);
  walker2.R[0] = TinyVector<double, 3>(1.0, 2.0, 0.0);

  crowd.addWalker(walker1);
  crowd.addWalker(walker2);

  std::vector<std::pair<int, int>> particle_group_indexes{{0, 1}, {1, 2}};

  RandomGenerator_t random_gen;

  MoveContext move_context(crowd.size(), num_particles, particle_group_indexes, random_gen);

  move_context.loadCrowd(crowd);

  using WalkerParticlePositions = std::vector<std::unique_ptr<MoveContext::ParticlePositions>>;
  auto returnWalkerPosition     = [](WalkerParticlePositions& vec_up_ppos, int w_index,
                                 int pos_index) -> TinyVector<double, 3> {
    return (*(vec_up_ppos[w_index]))[pos_index];
  };
  WalkerParticlePositions& walker_particle_positions = move_context.get_walker_positions();
  REQUIRE(returnWalkerPosition(walker_particle_positions, 0, 0) == TinyVector<double, 3>(1.0, 0.0, 0.0));
  REQUIRE(returnWalkerPosition(walker_particle_positions, 1, 0) == TinyVector<double, 3>(1.0, 2.0, 0.0));
}

} // namespace qmcplusplus
