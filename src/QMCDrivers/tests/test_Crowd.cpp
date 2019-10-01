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
#include "Estimators/tests/FakeEstimator.h"

#include "QMCApp/tests/MinimalWaveFunctionPool.h"
#include "QMCApp/tests/MinimalParticlePool.h"
#include "QMCApp/tests/MinimalHamiltonianPool.h"

namespace qmcplusplus
{
TEST_CASE("Crowd integration", "[drivers]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate* comm = OHMMS::Controller;

  EstimatorManagerBase em(comm);

  FakeEstimator* fake_est = new FakeEstimator;

  em.add(fake_est, "fake");

  ScalarEstimatorBase* est2 = em.getEstimator("fake");
  FakeEstimator* fake_est2  = dynamic_cast<FakeEstimator*>(est2);
  REQUIRE(fake_est2 != NULL);
  REQUIRE(fake_est2 == fake_est);

  Crowd crowd(em);

  
}
  
TEST_CASE("Crowd::loadWalkers", "[particle]")
{
  using WalkerMCP = Walker<QMCTraits, PtclOnLatticeTraits>;

  OHMMS::Controller->initialize(0, NULL);
  Communicate* comm = OHMMS::Controller;

  EstimatorManagerBase em(comm);
  FakeEstimator* fake_est = new FakeEstimator;
  em.add(fake_est, "fake");

  Crowd crowd(em);

  int num_particles = 1;
  WalkerMCP walker1(num_particles);
  TinyVector<double, 3> pos1(1.0, 0.0, 0.0);
  walker1.R[0] = pos1;
  WalkerMCP walker2(num_particles);
  TinyVector<double, 3> pos2(1.0, 2.0, 0.0);
  walker2.R[0] = pos2;

  MinimalParticlePool mpp;
  ParticleSetPool particle_pool = mpp(comm);
  MinimalWaveFunctionPool wfp(comm);
  WaveFunctionPool wavefunction_pool = wfp(&particle_pool);
  wavefunction_pool.setPrimary(wavefunction_pool.getWaveFunction("psi0"));
  MinimalHamiltonianPool mhp(comm);
  HamiltonianPool hamiltonian_pool = mhp(&particle_pool, &wavefunction_pool);
  
  ParticleSet pset1(*(particle_pool.getParticleSet("e")));
  pset1.R = pos1;
  ParticleSet pset2(*(particle_pool.getParticleSet("e")));
  pset2.R = pos2;

  std::unique_ptr<TrialWaveFunction> twf1;
  twf1.reset(wavefunction_pool.getPrimary()->makeClone(pset1));
  std::unique_ptr<TrialWaveFunction> twf2;
  twf2.reset(wavefunction_pool.getPrimary()->makeClone(pset2));
  std::unique_ptr<QMCHamiltonian> ham1;
  ham1.reset(hamiltonian_pool.getPrimary()->makeClone(pset1, *twf1));
  std::unique_ptr<QMCHamiltonian> ham2;
  ham1.reset(hamiltonian_pool.getPrimary()->makeClone(pset2, *twf2));
  
  crowd.addWalker(walker1,pset1,*twf1,*ham1);
  crowd.addWalker(walker2,pset2,*twf2,*ham2);

  std::vector<std::pair<int, int>> particle_group_indexes{{0, 1}, {1, 2}};

  RandomGenerator_t random_gen;

  crowd.loadWalkers();
  //actuall test loadWalkers  
}

TEST_CASE("Crowd::get_accept_ratio","[Drivers]")
{
  using WalkerMCP = Walker<QMCTraits, PtclOnLatticeTraits>;

  OHMMS::Controller->initialize(0, NULL);
  Communicate* comm = OHMMS::Controller;

  EstimatorManagerBase em(comm);
  FakeEstimator* fake_est = new FakeEstimator;
  em.add(fake_est, "fake");

  Crowd crowd(em);

  crowd.incAccept();
  crowd.incAccept();
  crowd.incAccept();
  crowd.incReject();
  REQUIRE( crowd.get_accept_ratio() == Approx(0.75));
}

}
