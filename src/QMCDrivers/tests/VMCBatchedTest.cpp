//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
#include "VMCBatchedTest.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "Concurrency/Info.hpp"
#include "CountingWaveFunctionComponent.h"
#include "EstimatorInputDelegates.h"
#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Particle/SampleStack.h"
#include "QMCDrivers/Crowd.h"
#include "QMCDrivers/MCPopulation.h"
#include "QMCDrivers/VMC/VMCBatched.h"
#include "QMCDrivers/VMC/VMCDriverInput.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/StdRandom.h"
#include <MinimalHamiltonianPool.h>
#include <MinimalParticlePool.h>
#include <MinimalWaveFunctionPool.h>

#include <cmath>
#include <limits>
#include <sstream>
#include <memory>
#include <string>
#include <vector>

namespace qmcplusplus::testing
{
VMCBatchedSummary VMCBatchedTest::runAllParticleNoDrift(int walkers_per_crowd)
{
  const int crowds           = Concurrency::maxCapacity<>();
  const int walkers_per_rank = crowds * walkers_per_crowd;
  const std::string input    = "<qmc method=\"vmc_batch\" move=\"allp\" checkpoint=\"-1\">"
                               "<parameter name=\"crowds\">" +
      std::to_string(crowds) +
      "</parameter>"
      "<parameter name=\"walkers_per_rank\">" +
      std::to_string(walkers_per_rank) +
      "</parameter>"
      "<parameter name=\"warmupsteps\">0</parameter><parameter name=\"substeps\">" +
      std::to_string(VMCBatchedTest::substeps) +
      "</parameter><parameter name=\"steps\">" +
      std::to_string(VMCBatchedTest::steps) +
      "</parameter><parameter name=\"blocks\">1</parameter>"
      "<parameter name=\"timestep\">0.01</parameter><parameter name=\"usedrift\">no</parameter></qmc>";
  ProjectData project_data("test", ProjectData::DriverVersion::BATCH);
  Communicate& comm = *OHMMS::Controller;
  Libxml2Document doc;
  REQUIRE(doc.parseFromString(input));
  xmlNodePtr node = doc.getRoot();
  QMCDriverInput qmcdriver_input;
  qmcdriver_input.readXML(node);
  VMCDriverInput vmcdriver_input;
  vmcdriver_input.readXML(node);
  // The pool randomizes initial electron positions with the global RNG, before any crowd tasks run.
  Random.init(7919 + comm.rank());
  auto particle_pool = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool =
      MinimalWaveFunctionPool::make_diamondC_1x1x1(project_data.getRuntimeOptions(), comm, particle_pool);
  wavefunction_pool.getWaveFunction().value().get().addComponent(std::make_unique<CountingWaveFunctionComponent>());
  auto hamiltonian_pool = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);

  using RNG = RandomBase<QMCTraits::FullPrecRealType>;
  std::vector<std::unique_ptr<StdRandom<QMCTraits::FullPrecRealType>>> rngs;
  RefVector<RNG> rng_refs;
  for (int crowd = 0; crowd < crowds; ++crowd)
  {
    const auto seed =
        static_cast<StdRandom<QMCTraits::FullPrecRealType>::uint_type>(7919 + 104729 * comm.rank() + 15485863 * crowd);
    rngs.push_back(std::make_unique<StdRandom<QMCTraits::FullPrecRealType>>(seed));
    rng_refs.push_back(*rngs.back());
  }
  WalkerConfigurations walker_confs;
  SampleStack samples;
  VMCBatched driver(project_data, std::move(qmcdriver_input), nullptr, std::move(vmcdriver_input), walker_confs,
                    MCPopulation(comm.size(), comm.rank(), *particle_pool.getParticleSet("e"),
                                 wavefunction_pool.getWaveFunction().value(),
                                 hamiltonian_pool.getHamiltonian().value()),
                    rng_refs, samples, comm);
  driver.setStatus("vmcbatched_allp", "", false);
  driver.process(node);
  driver.run();

  VMCBatchedSummary summary;
  summary.walkers = driver.get_num_living_walkers();
  const auto fail = [&summary](const std::string& category, int crowd, int walker = -1) {
    std::ostringstream message;
    message << category;
    if (crowd >= 0)
    {
      message << " (crowd " << crowd;
      if (walker >= 0)
        message << ", walker " << walker;
      message << ')';
    }
    summary.failures.push_back(message.str());
  };
  if (driver.crowds_.size() != static_cast<std::size_t>(crowds) || summary.walkers != walkers_per_rank)
    fail("distribution mismatch", -1);
  for (std::size_t crowd_index = 0; crowd_index < driver.crowds_.size(); ++crowd_index)
  {
    const auto& crowd = driver.crowds_[crowd_index];
    if (crowd->size() != static_cast<std::size_t>(walkers_per_crowd))
      fail("distribution mismatch", crowd_index);
    summary.particle_moves += crowd->get_accept() + crowd->get_reject();
    const auto& walker_twfs  = crowd->get_walker_twfs();
    const auto& walker_elecs = crowd->get_walker_elecs();
    for (const auto& particles : walker_elecs)
      summary.expected_particle_moves += VMCBatchedTest::moves_per_walker * particles.get().getTotalNum();
    if (walker_twfs.size() != walker_elecs.size())
      fail("missing component", crowd_index);
    for (std::size_t walker = 0; walker < walker_twfs.size(); ++walker)
    {
      if (walker >= walker_elecs.size())
      {
        fail("missing component", crowd_index, walker);
        continue;
      }
      const auto& wf       = walker_twfs[walker];
      const auto& orbitals = wf.get().getOrbitals();
      const auto* counter =
          orbitals.empty() ? nullptr : dynamic_cast<const CountingWaveFunctionComponent*>(orbitals.back().get());
      if (!counter || !std::isfinite(wf.get().getLogPsi()))
      {
        fail("missing component or nonfinite log", crowd_index, walker);
        continue;
      }

      const auto& proposals   = counter->getProposals();
      const auto& resolutions = counter->getResolutions();
      summary.proposals += proposals.size();
      summary.resolutions += resolutions.size();
      if (proposals.size() != static_cast<std::size_t>(VMCBatchedTest::moves_per_walker) ||
          resolutions.size() != static_cast<std::size_t>(VMCBatchedTest::moves_per_walker))
        fail("record count mismatch", crowd_index, walker);
      for (const auto& proposal : proposals)
        if (!proposal.evaluate || proposal.evaluations != 1)
          fail("proposal mask/evaluation mismatch", crowd_index, walker);
      for (const auto& resolution : resolutions)
      {
        if (resolution.accepted)
          ++summary.accepted;
        else
          ++summary.rejected;
        if (resolution.mixed)
          ++summary.mixed_resolutions;
        if (resolution.evaluations != static_cast<int>(!resolution.accepted))
          fail("resolution record/evaluation mismatch", crowd_index, walker);
      }

      // Check the cached final state against a fresh serial evaluation at unchanged positions.
      ParticleSet& particles = walker_elecs[walker].get();
      const auto saved_log   = wf.get().getLogPsi();
      const auto saved_G     = particles.G;
      const auto saved_L     = particles.L;
      const auto fresh_log   = wf.get().evaluateLog(particles);
      // Use the working precision of the wavefunction calculation, including mixed-precision builds.
      const auto tolerance = std::numeric_limits<QMCTraits::RealType>::epsilon() * 1000;
      const auto close     = [tolerance](auto actual, auto expected) {
        return std::abs(actual - expected) <= tolerance * (1 + std::abs(expected));
      };
      if (!close(saved_log, fresh_log) || saved_G.size() != particles.G.size() || saved_L.size() != particles.L.size())
        fail("log/GL numeric mismatch", crowd_index, walker);
      else
      {
        for (std::size_t i = 0; i < saved_G.size(); ++i)
          for (int d = 0; d < OHMMS_DIM; ++d)
            if (!close(saved_G[i][d], particles.G[i][d]))
              fail("GL numeric mismatch", crowd_index, walker);
        for (std::size_t i = 0; i < saved_L.size(); ++i)
          if (!close(saved_L[i], particles.L[i]))
            fail("L numeric mismatch", crowd_index, walker);
      }
    }
  }
  if (summary.particle_moves != summary.expected_particle_moves)
    fail("particle move count mismatch", -1);
  if (summary.accepted + summary.rejected != summary.resolutions)
    fail("resolution accounting mismatch", -1);
  return summary;
}

} // namespace qmcplusplus::testing
