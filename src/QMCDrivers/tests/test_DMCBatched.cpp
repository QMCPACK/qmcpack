//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
#include <string>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include "Utilities/for_testing/Catch2Approx.h"

#include "Message/Communicate.h"
#include "Message/UniformCommunicateError.h"
#include "QMCDrivers/DMC/DMCDriverInput.h"
#include "QMCDrivers/DMC/DMCBatched.h"
#include "QMCDrivers/DMC/WalkerControl.h"
#include "QMCDrivers/SFNBranch.h"
#include "QMCDrivers/Crowd.h"
#include "Estimators/EstimatorManagerNew.h"
#include "Particle/HDFWalkerOutput.h"
#include "Particle/HDFWalkerInputManager.h"
#include "QMCDrivers/QMCDriverFactory.h"
#include "QMCDrivers/tests/ValidQMCInputSections.h"
#include "QMCDrivers/tests/SetupDMCTest.h"
#include "EstimatorInputDelegates.h"
#include "Concurrency/Info.hpp"
#include "Concurrency/UtilityFunctions.hpp"
#include "Platforms/Host/OutputManager.h"
#include "SetupPools.h"

namespace qmcplusplus
{
namespace testing
{
class DMCBatchedTest
{
public:
  DMCBatchedTest() { up_dtest_ = std::make_unique<SetupDMCTest>(1); }

  /** Save live population weights through the production configuration/checkpoint path. */
  static void checkpoint(DMCBatched& driver, const std::string& root, int block)
  {
    driver.population_.saveWalkerConfigurations(driver.walker_configs_ref_);
    driver.setWalkerOffsets(driver.walker_configs_ref_, driver.myComm);
    HDFWalkerOutput output(driver.population_.get_golden_electrons().getTotalNum(), driver.myComm);
    output.dump(driver.walker_configs_ref_, root, block);
  }

  /** Verify initialization and actual DMC propagation, skipped/conventional combs, and block boundaries. */
  static void checkWeightLifecycle(DMCBatched& driver)
  {
    const bool fixed = driver.dmcdriver_input_.get_reconfiguration();
    auto& population = driver.population_;
    const int count  = population.get_num_global_walkers();
    std::vector<double> expected_weights;
    for (const auto& walker : population.get_walkers())
      expected_weights.push_back(fixed ? walker->Weight : 1);
    auto contexts = driver.getContextForStepsRefs();
    for (int crowd = 0; crowd < driver.crowds_.size(); ++crowd)
      driver.initialLogEvaluation(crowd, driver.crowds_, contexts, driver.serializing_crowd_walkers_, fixed);
    for (int iw = 0; iw < expected_weights.size(); ++iw)
      CHECK(population.get_walkers()[iw]->Weight == expected_weights[iw]);

    double energy, variance;
    population.measureGlobalEnergyVariance(driver.myComm, energy, variance, fixed);
    driver.branch_engine_->initParam(population, energy, variance, fixed, false);
    DMCBatched::StateForThread state(driver.qmcdriver_input_, *driver.drift_modifier_, *driver.branch_engine_,
                                     population, 2, driver.serializing_crowd_walkers_);
    state.recalculate_properties_period = 1;
    driver.estimator_manager_->startDriverRun(driver.get_root_name());
    for (int block = 0; block < 2; ++block)
    {
      std::vector<double> measured_moments(3, 0);
      expected_weights.clear();
      for (const auto& walker : population.get_walkers())
        expected_weights.push_back(walker->Weight);
      driver.estimator_manager_->startBlock(2);
      for (auto& crowd : driver.crowds_)
        crowd->startBlock(2);
      for (int iw = 0; iw < expected_weights.size(); ++iw)
        CHECK(population.get_walkers()[iw]->Weight == expected_weights[iw]);
      for (int step = 0; step < 2; ++step)
      {
        std::vector<double> old_weights, old_energies;
        for (const auto& walker : population.get_walkers())
        {
          old_weights.push_back(walker->Weight);
          old_energies.push_back(walker->Properties(WalkerProperties::Indexes::LOCALENERGY));
        }
        state.step        = step;
        state.global_step = block * 2 + step;
        for (int crowd = 0; crowd < driver.crowds_.size(); ++crowd)
          driver.runDMCStep(crowd, state, driver.timers_, driver.dmc_timers_, driver.step_contexts_, driver.crowds_);
        double weight = 0;
        std::vector<double> propagated_weights;
        for (int iw = 0; iw < population.get_num_local_walkers(); ++iw)
        {
          const auto& walker = *population.get_walkers()[iw];
          const auto factor =
              driver.branch_engine_->branchWeight(walker.Properties(WalkerProperties::Indexes::LOCALENERGY),
                                                  old_energies[iw]);
          CHECK(walker.Weight == Approx(old_weights[iw] * factor));
          weight += walker.Weight;
          propagated_weights.push_back(walker.Weight);
          if (step % driver.qmcdriver_input_.get_estimator_measurement_period() == 0)
          {
            const double energy = walker.Properties(WalkerProperties::Indexes::LOCALENERGY);
            measured_moments[0] += walker.Weight;
            measured_moments[1] += walker.Weight * energy;
            measured_moments[2] += walker.Weight * energy * energy;
          }
        }
        driver.myComm.allreduce(weight);
        const bool skip = state.global_step == 0;
        driver.walker_controller_->branch(state.global_step, population, skip);
        if (fixed)
        {
          CHECK(population.get_num_global_walkers() == count);
          CHECK(driver.walker_controller_->get_ensemble_property().Weight == Approx(weight));
          if (skip)
            for (int iw = 0; iw < propagated_weights.size(); ++iw)
              CHECK(population.get_walkers()[iw]->Weight == propagated_weights[iw]);
          else
            for (const auto& walker : population.get_walkers())
              CHECK(walker->Weight == Approx(weight / count));
        }
        driver.branch_engine_->updateParamAfterPopControl(driver.walker_controller_->get_ensemble_property(),
                                                          population.get_golden_electrons().getTotalNum());
        population.redistributeWalkers(driver.crowds_);
      }
      driver.endBlock();
      driver.myComm.allreduce(measured_moments);
      if (driver.myComm.rank() == 0)
      {
        CHECK(driver.estimator_manager_->get_AverageCache()[0] == Approx(measured_moments[1] / measured_moments[0]));
        CHECK(driver.estimator_manager_->get_AverageCache()[1] == Approx(measured_moments[2] / measured_moments[0]));
      }
    }
    driver.estimator_manager_->stopDriverRun();
    population.saveWalkerConfigurations(driver.walker_configs_ref_);
    // The shared VMC/dynamic-DMC initialization contract remains an explicit reset to unit weights.
    for (int crowd = 0; crowd < driver.crowds_.size(); ++crowd)
      driver.initialLogEvaluation(crowd, driver.crowds_, contexts, driver.serializing_crowd_walkers_, false);
    for (const auto& walker : population.get_walkers())
      CHECK(walker->Weight == 1);
  }

private:
  UPtr<SetupDMCTest> up_dtest_;
};
} // namespace testing

TEST_CASE("DMCDriverInput L2 diffusion", "[drivers]")
{
  DMCDriverInput dmcdriver_input;
  CHECK_FALSE(dmcdriver_input.get_l2_diffusion());

  Libxml2Document doc;
  REQUIRE(doc.parseFromString(R"(<qmc method="dmc"><parameter name="L2_diffusion">yes</parameter></qmc>)"));
  dmcdriver_input.readXML(doc.getRoot());
  CHECK(dmcdriver_input.get_l2_diffusion());
}

TEST_CASE("DMCDriverInput stochastic comb opt-in and compatibility", "[drivers][stochastic_comb]")
{
  CHECK_FALSE(DMCDriverInput{}.get_reconfiguration());
  for (const std::string option : {"", "no", "yes", "runwhileincorrect"})
  {
    CAPTURE(option);
    DMCDriverInput input;
    Libxml2Document doc;
    const std::string xml =
        option.empty() ? "<qmc/>" : "<qmc><parameter name=\"reconfiguration\">" + option + "</parameter></qmc>";
    REQUIRE(doc.parseFromString(xml));
    input.readXML(doc.getRoot());
    CHECK(input.get_reconfiguration() == (option == "yes" || option == "runwhileincorrect"));
  }
}

TEST_CASE("DMCBatched fixed weight initialization restart and real propagation", "[drivers][stochastic_comb]")
{
  using namespace testing;
  const auto reconfiguration = GENERATE("yes", "runwhileincorrect");
  CAPTURE(reconfiguration);
  Communicate& comm(*OHMMS::Controller);
  ProjectData project;
  SetupPools pools;
  RandomNumberGeneratorPool rngs(2);
  WalkerConfigurations configs;
  auto& electrons = *pools.particle_pool->getParticleSet("e");
  configs.createWalkers(4, electrons.getTotalNum());
  for (int iw = 0; iw < 4; ++iw)
  {
    configs[iw]->R = electrons.R;
    configs[iw]->R[0][0] += 0.02 * (iw + 4 * comm.rank());
    configs[iw]->Weight = 0.5 * (iw + 1);
  }
  Libxml2Document doc;
  REQUIRE(doc.parseFromString(std::string(R"(<qmc method="dmc" move="pbyp">
    <parameter name="reconfiguration">)") +
                              reconfiguration + R"(</parameter>
    <parameter name="walkers_per_rank">4</parameter>
    <parameter name="crowds">2</parameter>
    <parameter name="warmupsteps">1</parameter>
    <parameter name="blocks">2</parameter>
    <parameter name="steps">2</parameter>
    <parameter name="timestep">0.01</parameter>
    <parameter name="target_walkers">24</parameter>
    <parameter name="estimator_period">2</parameter>
  </qmc>)"));
  // Construct a new DMC section from either in-memory or restored configurations.
  auto makeDriver = [&](WalkerConfigurations& input, const std::string& root) {
    QMCDriverInput qmc_input;
    qmc_input.readXML(doc.getRoot());
    DMCDriverInput dmc_input;
    dmc_input.readXML(doc.getRoot());
    auto driver = std::make_unique<DMCBatched>(project, std::move(qmc_input), nullptr, std::move(dmc_input), input,
                                               MCPopulation(comm.size(), comm.rank(), electrons,
                                                            pools.wavefunction_pool->getWaveFunction().value(),
                                                            pools.hamiltonian_pool->getHamiltonian().value()),
                                               rngs.getRngRefs(), comm);
    driver->setStatus(root, "", false);
    driver->process(doc.getRoot());
    return driver;
  };
  auto driver                       = makeDriver(configs, "fixed_dmc_cpp_original");
  const std::string checkpoint_root = "fixed_dmc_cpp_restart";
  DMCBatchedTest::checkpoint(*driver, checkpoint_root, 1);
  comm.barrier();
  WalkerConfigurations restored;
  HDFWalkerInputManager reader(restored, electrons.getTotalNum(), comm);
  Libxml2Document input_doc;
  REQUIRE(input_doc.parseFromString("<mcwalkerset fileroot=\"" + checkpoint_root + "\" version=\"0 4\"/>"));
  REQUIRE(reader.put(input_doc.getRoot()));
  REQUIRE(restored.getActiveWalkers() == configs.getActiveWalkers());
  for (int iw = 0; iw < restored.getActiveWalkers(); ++iw)
  {
    CHECK(restored[iw]->Weight == configs[iw]->Weight);
    CHECK(restored[iw]->R[0][0] == configs[iw]->R[0][0]);
  }
  auto restarted = makeDriver(restored, "fixed_dmc_cpp_restored");
  std::vector<RandomBase<double>::uint_type> master_state;
  Random.save(master_state);
  DMCBatchedTest::checkWeightLifecycle(*driver);
  Random.load(master_state);
  DMCBatchedTest::checkWeightLifecycle(*restarted);
  for (int iw = 0; iw < restored.getActiveWalkers(); ++iw)
  {
    CHECK(restored[iw]->Weight == Approx(configs[iw]->Weight));
    CHECK(restored[iw]->R[0][0] == configs[iw]->R[0][0]);
  }
  auto next_section = makeDriver(configs, "fixed_dmc_cpp_next_section");
  DMCBatchedTest::checkWeightLifecycle(*next_section);
}

TEST_CASE("QMCDriverFactory rejects invalid L2 diffusion configurations for DMCBatched", "[drivers]")
{
  using namespace testing;
  RandomNumberGeneratorPool rng_pool(1);
  ProjectData test_project;
  Communicate& comm(*OHMMS::Controller);

  Libxml2Document doc;
  REQUIRE(doc.parseFromString(R"(<qmc method="dmc"><parameter name="L2_diffusion">yes</parameter></qmc>)"));
  xmlNodePtr node = doc.getRoot();

  auto particle_pool = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool =
      MinimalWaveFunctionPool::make_diamondC_1x1x1(test_project.getRuntimeOptions(), comm, particle_pool);
  auto hamiltonian_pool = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);

  std::string expected_error;
  SECTION("spinor ParticleSet")
  {
    particle_pool.getParticleSet("e")->setSpinor(true);
    expected_error = "L2 diffusion is not supported for spinor particle sets.";
  }
  SECTION("Hamiltonian without L2 potential")
  { expected_error = "L2 diffusion was requested, but the Hamiltonian has no L2 potential."; }

  QMCDriverFactory driver_factory(test_project);
  auto das = driver_factory.readSection(node);

  auto construct_driver = [&]() {
    driver_factory.createQMCDriver(node, das, std::nullopt, *particle_pool.getWalkerSet("e"), particle_pool,
                                   wavefunction_pool, hamiltonian_pool, comm);
  };

  CHECK_THROWS_MATCHES(construct_driver(), UniformCommunicateError, Catch::Matchers::Message(expected_error));
}

/** Since we check the DMC only feature of reserve walkers perhaps this should be
 *  a DMC integration test.
 */
#ifdef _OPENMP
TEST_CASE("DMCDriver+QMCDriverNew integration", "[drivers]")
{
  using namespace testing;
  Concurrency::OverrideMaxCapacity<> override(8);
  RandomNumberGeneratorPool rng_pool(8);
  ProjectData test_project;
  Communicate& comm(*OHMMS::Controller);
  outputManager.pause();

  Libxml2Document doc;
  REQUIRE(doc.parseFromString(valid_dmc_input_sections[valid_dmc_input_dmc_batch_index]));
  xmlNodePtr node = doc.getRoot();
  QMCDriverInput qmcdriver_input;
  qmcdriver_input.readXML(node);
  DMCDriverInput dmcdriver_input;
  dmcdriver_input.readXML(node);
  auto particle_pool = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool =
      MinimalWaveFunctionPool::make_diamondC_1x1x1(test_project.getRuntimeOptions(), comm, particle_pool);

  auto hamiltonian_pool = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);
  SampleStack samples;
  WalkerConfigurations walker_confs;

  DMCBatched dmcdriver(test_project, std::move(qmcdriver_input), nullptr, std::move(dmcdriver_input), walker_confs,
                       MCPopulation(comm.size(), comm.rank(), *particle_pool.getParticleSet("e"),
                                    wavefunction_pool.getWaveFunction().value(),
                                    hamiltonian_pool.getHamiltonian().value()),
                       rng_pool.getRngRefs(), comm);

  // setStatus must be called before process
  std::string root_name{"Test"};
  //For later sections this appears to contain important state.
  std::string prev_config_file{""};

  dmcdriver.setStatus(root_name, prev_config_file, false);
  // We want to express out expectations of the QMCDriver state machine so we catch
  // changes to it over time.
  outputManager.resume();

  dmcdriver.process(node);
  CHECK(dmcdriver.get_num_living_walkers() == 8 / comm.size());
  const QMCTraits::IndexType reserved_walkers = dmcdriver.get_num_living_walkers() + dmcdriver.get_num_dead_walkers();
  CHECK(reserved_walkers == 10 / comm.size());
  // What else should we expect after process
}
#endif

} // namespace qmcplusplus
