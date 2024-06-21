//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Message/Communicate.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Particle/tests/MinimalParticlePool.h"
#include "QMCHamiltonians/tests/MinimalHamiltonianPool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
#include "Utilities/ResourceCollection.h"
#include "Utilities/ProjectData.h"

#include "QMCDrivers/WalkerLogCollector.h"


namespace qmcplusplus
{

TEST_CASE("WalkerLogCollector::collect", "[estimators]")
{
  app_log() << "\n\n=======================================================\n";
  app_log() << "test WalkerLogCollector::collect\n";

  ProjectData test_project("test", ProjectData::DriverVersion::BATCH);
  Communicate* comm = OHMMS::Controller;


  auto particle_pool = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool =
      MinimalWaveFunctionPool::make_diamondC_1x1x1(test_project.getRuntimeOptions(), comm, particle_pool);
  auto& pset = *(particle_pool.getParticleSet("e"));
  // This is where the pset properties "properies" gain the different hamiltonian operator values.
  auto hamiltonian_pool = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);

  auto& twf = *(wavefunction_pool.getWaveFunction("wavefunction"));
  auto& ham = *(hamiltonian_pool.getPrimary());

  // setup data structures for multiple walkers

  UPtrVector<QMCHamiltonian> hams;
  UPtrVector<TrialWaveFunction> twfs;
  std::vector<ParticleSet> psets;

  int num_walkers   = 4;
  int num_electrons = particle_pool.getParticleSet("e")->getTotalNum();
  int num_ions      = particle_pool.getParticleSet("ion")->getTotalNum();

  for (int iw = 0; iw < num_walkers; ++iw)
  {
    psets.emplace_back(pset);
    psets.back().randomizeFromSource(*particle_pool.getParticleSet("ion"));
    twfs.emplace_back(twf.makeClone(psets.back()));
    hams.emplace_back(hamiltonian_pool.getPrimary()->makeClone(psets.back(), *twfs.back()));
  }

  using MCPWalker = Walker<QMCTraits, PtclOnLatticeTraits>;

  std::vector<MCPWalker> walkers(num_walkers, MCPWalker(pset.getTotalNum()));

  for (auto& walker : walkers)
  {
    walker.R          = pset.R;
    walker.spins      = pset.spins;
    walker.Properties = pset.Properties;
    walker.registerData();
    walker.DataSet.allocate();
  }

  WalkerLogState state{true, 1, false};
  WalkerLogCollector wlc(state);

  auto& bsi = wlc.walker_property_int_buffer;
  auto& bsr = wlc.walker_property_real_buffer;
  auto& bpr = wlc.walker_particle_real_buffer;

  CHECK(bsi.nrows() == 0);
  CHECK(bsr.nrows() == 0);
  CHECK(bpr.nrows() == 0);
  CHECK(bsi.ncols() == 0);
  CHECK(bsr.ncols() == 0);
  CHECK(bpr.ncols() == 0);

#ifndef QMC_COMPLEX
  const size_t npcols = 56;
#else
  const size_t npcols = 88;
#endif

  int step = 0;
  wlc.collect(walkers[0], psets[0], *(twfs[0]), *(hams[0]), step);

  CHECK(bsi.nrows() == 1);
  CHECK(bsr.nrows() == 1);
  CHECK(bpr.nrows() == 1);
  CHECK(bsi.ncols() == 4);
  CHECK(bsr.ncols() == 13);
  CHECK(bpr.ncols() == npcols);

  for (size_t iw = 1; iw < walkers.size(); ++iw)
    wlc.collect(walkers[iw], psets[iw], *(twfs[iw]), *(hams[iw]), step);

  CHECK(bsi.nrows() == 4);
  CHECK(bsr.nrows() == 4);
  CHECK(bpr.nrows() == 4);
  CHECK(bsi.ncols() == 4);
  CHECK(bsr.ncols() == 13);
  CHECK(bpr.ncols() == npcols);

  for (step = 1; step < 3; ++step)
    for (size_t iw = 0; iw < walkers.size(); ++iw)
      wlc.collect(walkers[iw], psets[iw], *(twfs[iw]), *(hams[iw]), step);

  CHECK(bsi.nrows() == 12);
  CHECK(bsr.nrows() == 12);
  CHECK(bpr.nrows() == 12);
  CHECK(bsi.ncols() == 4);
  CHECK(bsr.ncols() == 13);
  CHECK(bpr.ncols() == npcols);

  app_log() << "\nend test WalkerLogCollector::collect\n";
  app_log() << "=======================================================\n";
}


} // namespace qmcplusplus
