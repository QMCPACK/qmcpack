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
#include "Message/Communicate.h"
#include "QMCDrivers/Crowd.h"
#include "type_traits/template_types.hpp"
#include "Estimators/tests/FakeEstimator.h"
#include "Estimators/EstimatorManagerNew.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
#include "Particle/tests/MinimalParticlePool.h"
#include "QMCHamiltonians/tests/MinimalHamiltonianPool.h"

#include "QMCDrivers/tests/SetupPools.h"

namespace qmcplusplus
{
namespace testing
{
class CrowdWithWalkers
{
public:
  using MCPWalker = Walker<QMCTraits, PtclOnLatticeTraits>;

  EstimatorManagerNew em;
  UPtr<Crowd> crowd_ptr;
  Crowd& get_crowd() { return *crowd_ptr; }
  UPtrVector<MCPWalker> walkers;
  UPtrVector<ParticleSet> psets;
  UPtrVector<TrialWaveFunction> twfs;
  UPtrVector<QMCHamiltonian> hams;
  std::vector<TinyVector<double, 3>> tpos;
  DriverWalkerResourceCollection fatwalker_resource_collection_;
  const MultiWalkerDispatchers dispatchers_;

public:
  CrowdWithWalkers(SetupPools& pools) : em(pools.comm), dispatchers_(true)
  {
    FakeEstimator* fake_est = new FakeEstimator;
    em.add(fake_est, "fake");

    crowd_ptr    = std::make_unique<Crowd>(em, fatwalker_resource_collection_, dispatchers_);
    Crowd& crowd = *crowd_ptr;
    // To match the minimal particle set
    int num_particles = 2;
    // for testing we update the first position in the walker
    auto makePointWalker = [this, &pools, &crowd, num_particles](TinyVector<double, 3> pos) {
      walkers.emplace_back(std::make_unique<MCPWalker>(num_particles));
      walkers.back()->R[0] = pos;
      psets.emplace_back(std::make_unique<ParticleSet>(*(pools.particle_pool->getParticleSet("e"))));
      twfs.emplace_back(pools.wavefunction_pool->getPrimary()->makeClone(*psets.back()));
      hams.emplace_back(pools.hamiltonian_pool->getPrimary()->makeClone(*psets.back(), *twfs.back()));
      crowd.addWalker(*walkers.back(), *psets.back(), *twfs.back(), *hams.back());
    };

    tpos.push_back(TinyVector<double, 3>(1.0, 0.0, 0.0));
    makePointWalker(tpos.back());
    tpos.push_back(TinyVector<double, 3>(1.0, 2.0, 0.0));
    makePointWalker(tpos.back());
  }

  void makeAnotherPointWalker()
  {
    walkers.emplace_back(std::make_unique<MCPWalker>(*walkers.back()));
    psets.emplace_back(std::make_unique<ParticleSet>(*psets.back()));
    twfs.emplace_back(twfs.back()->makeClone(*psets.back()));
    hams.emplace_back(hams.back()->makeClone(*psets.back(), *twfs.back()));
  }
};
} // namespace testing

TEST_CASE("Crowd integration", "[drivers]")
{
  Communicate* comm = OHMMS::Controller;

  EstimatorManagerNew em(comm);

  FakeEstimator* fake_est = new FakeEstimator;

  em.add(fake_est, "fake");

  ScalarEstimatorBase* est2 = em.getEstimator("fake");
  FakeEstimator* fake_est2  = dynamic_cast<FakeEstimator*>(est2);
  REQUIRE(fake_est2 != nullptr);
  REQUIRE(fake_est2 == fake_est);
  // The above was required behavior for crowd at one point.
  // TODO: determine whether it still is, I don't think so.
  const MultiWalkerDispatchers dispatchers(true);
  DriverWalkerResourceCollection fatwalker_resource_collection_;
  Crowd crowd(em, fatwalker_resource_collection_, dispatchers);
}

TEST_CASE("Crowd::loadWalkers", "[particle]")
{
  using namespace testing;
  SetupPools pools;

  CrowdWithWalkers crowd_with_walkers(pools);
  Crowd& crowd = crowd_with_walkers.get_crowd();

  std::vector<std::pair<int, int>> particle_group_indexes{{0, 1}, {1, 2}};

  RandomGenerator_t random_gen;

  crowd.loadWalkers();

  auto checkParticleSetPos = [&crowd_with_walkers](int iw) {
    REQUIRE(crowd_with_walkers.psets[iw]->R[0] == crowd_with_walkers.tpos[iw]);
  };
  for (int i = 0; i < crowd.size(); ++i)
    checkParticleSetPos(i);
}

TEST_CASE("Crowd redistribute walkers")
{
  using namespace testing;
  SetupPools pools;

  CrowdWithWalkers crowd_with_walkers(pools);
  Crowd& crowd = crowd_with_walkers.get_crowd();

  crowd_with_walkers.makeAnotherPointWalker();
  crowd.clearWalkers();
  for (int iw = 0; iw < crowd_with_walkers.walkers.size(); ++iw)
    crowd.addWalker(*crowd_with_walkers.walkers[iw], *crowd_with_walkers.psets[iw], *crowd_with_walkers.twfs[iw],
                    *crowd_with_walkers.hams[iw]);
  REQUIRE(crowd.size() == 3);
}

} // namespace qmcplusplus
