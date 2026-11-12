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
#include "Estimators/EstimatorManagerNew.h"
#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Particle/DistanceTable.h"
#include "Particle/ParticleBase/RandomSeqGenerator.h"
#include "Particle/SampleStack.h"
#include "QMCDrivers/Crowd.h"
#include "QMCDrivers/GreenFunctionModifiers/DriftModifierUNR.h"
#include "QMCDrivers/MCPopulation.h"
#include "QMCDrivers/VMC/VMCBatched.h"
#include "QMCDrivers/VMC/VMCDriverInput.h"
#include "QMCHamiltonians/BareKineticEnergy.h"
#include "Utilities/FakeRandom.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/StdRandom.h"
#include "Utilities/for_testing/Catch2Approx.h"
#include <MinimalHamiltonianPool.h>
#include <MinimalParticlePool.h>
#include <MinimalWaveFunctionPool.h>

#include <cmath>
#include <limits>
#include <memory>
#include <sstream>
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
  // Match QMCMain setup: initialize distance tables and structure factors before cloning walkers.
  particle_pool.reset();

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

void VMCBatchedTest::checkAllParticleBoundaryRejection(int invalid_particle)
{
  Libxml2Document doc;
  REQUIRE(doc.parseFromString(R"(<qmc method="vmc_batch" move="allp">
    <parameter name="substeps">1</parameter>
    <parameter name="timestep">0.01</parameter>
    <parameter name="usedrift">no</parameter></qmc>)"));
  QMCDriverInput qmc_input;
  qmc_input.readXML(doc.getRoot());
  VMCDriverInput vmc_input;
  vmc_input.readXML(doc.getRoot());

  Lattice lattice;
  lattice.BoxBConds          = false;
  lattice.R                  = ParticleSet::Tensor_t(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
  lattice.explicitly_defined = true;
  lattice.reset();
  const SimulationCell cell(lattice);
  ParticleSet elecs(cell);
  elecs.setName("e");
  elecs.create({2});
  elecs.R[0]              = {0.5, 0.5, 0.5};
  elecs.R[1]              = {0.7, 0.5, 0.5};
  auto& species           = elecs.getSpeciesSet();
  const int electron      = species.addSpecies("u");
  const int mass          = species.addAttribute("mass");
  species(mass, electron) = 1.0;
  elecs.resetGroups();
  const int aa_table = elecs.addTable(elecs);
  elecs.update();

  RuntimeOptions runtime_options;
  TrialWaveFunction psi(runtime_options);
  psi.addComponent(std::make_unique<CountingWaveFunctionComponent>());
  QMCHamiltonian ham;
  ham.addOperator(std::make_unique<BareKineticEnergy>(elecs), "Kinetic");
  ham.addObservables(elecs);
  Communicate& comm = *OHMMS::Controller;
  MCPopulation population(comm.size(), comm.rank(), elecs, psi, ham);
  WalkerConfigurations configs;
  population.createWalkers(2, configs);

  DriverWalkerResourceCollection resources;
  elecs.createResource(resources.pset_res);
  psi.createResource(resources.twf_res);
  ham.createResource(resources.ham_res);
  EstimatorManagerNew estimators(ham, comm);
  Crowd crowd(estimators, resources, elecs, psi, ham);

  // With FakeRandom(0.5), the particle-major displacements in walker 1
  // are (0, -sqrt(2*log(2)*tau), 0) for both particles. Only one crosses y=0.
  population.get_elec_particle_sets()[1]->R[invalid_particle][1] = 0.05;
  std::vector<ParticleSet::ParticlePos> initial_positions;
  for (int iw = 0; iw < 2; ++iw)
  {
    auto& particles = *population.get_elec_particle_sets()[iw];
    auto& twf       = *population.get_twfs()[iw];
    auto& walker    = *population.get_walkers()[iw];
    particles.update();
    twf.evaluateLog(particles);
    particles.saveWalker(walker);
    initial_positions.push_back(particles.R);
    crowd.addWalker(walker, particles, twf, *population.get_hamiltonians()[iw]);
  }

  FakeRandom<QMCTraits::FullPrecRealType> preview_rng;
  MCCoords<CoordsType::POS> deltas(4);
  makeGaussRandomWithEngine(deltas, preview_rng);
  for (int ip = 0; ip < 2; ++ip)
    for (int iw = 0; iw < 2; ++iw)
    {
      const int index = ip * 2 + iw;
      deltas.positions[index] *= std::sqrt(qmc_input.get_tau());
      const bool valid = lattice.isValid(lattice.toUnit(initial_positions[iw][ip] + deltas.positions[index]));
      REQUIRE(valid == !(iw == 1 && ip == invalid_particle));
    }

  DriftModifierUNR drift;
  VMCBatched::StateForThread state(qmc_input, vmc_input, drift, population, 1, false);
  QMCDriverNew::DriverTimers timers("VMCBatchedBoundaryTest::");
  FakeRandom<QMCTraits::FullPrecRealType> rng;
  QMCDriverNew::ContextForSteps context(rng);
  // The constant wavefunction gives probability one, so only boundary validity can reject.
  VMCBatched::advanceWalkers<CoordsType::POS>(state, crowd, timers, context, false, false);

  CHECK(crowd.get_accept() == 2);
  CHECK(crowd.get_reject() == 2);
  for (int iw = 0; iw < 2; ++iw)
  {
    CAPTURE(iw, invalid_particle);
    const auto& particles = *population.get_elec_particle_sets()[iw];
    const auto& walker    = *population.get_walkers()[iw];
    auto& twf             = *population.get_twfs()[iw];
    const auto& counter   = dynamic_cast<const CountingWaveFunctionComponent&>(*twf.getOrbitals().front());
    REQUIRE(counter.getProposals().size() == 1);
    CHECK(counter.getProposals()[0].evaluate);
    CHECK(counter.getProposals()[0].evaluations == 1);
    REQUIRE(counter.getResolutions().size() == 1);
    CHECK(counter.getResolutions()[0].accepted == (iw == 0));
    CHECK(counter.getResolutions()[0].evaluations == (iw == 1 ? 1 : 0));
    CHECK(counter.getResolutions()[0].mixed);
    CHECK(twf.getLogPsi() == Approx(0.0));
    auto expected = initial_positions[iw];
    for (int ip = 0; ip < 2; ++ip)
    {
      if (iw == 0)
        expected[ip] += deltas.positions[ip * 2 + iw];
      for (int d = 0; d < OHMMS_DIM; ++d)
      {
        CHECK(particles.R[ip][d] == Approx(expected[ip][d]));
        CHECK(particles.getCoordinates().getAllParticlePos()[ip][d] == Approx(expected[ip][d]));
        CHECK(walker.R[ip][d] == Approx(expected[ip][d]));
      }
    }
    const ParticleSet::PosType separation = expected[1] - expected[0];
    CHECK(particles.getDistTableAA(aa_table).getDistances()[1][0] == Approx(std::sqrt(dot(separation, separation))));
  }
}

} // namespace qmcplusplus::testing
