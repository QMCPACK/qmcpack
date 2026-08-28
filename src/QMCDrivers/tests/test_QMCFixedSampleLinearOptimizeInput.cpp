//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>

#include "EstimatorInputDelegates.h"
#include "QMCDrivers/WFOpt/QMCFixedSampleLinearOptimize.h"
#include "QMCDrivers/WFOpt/QMCFixedSampleLinearOptimizeBatched.h"
#include "QMCDrivers/MCPopulation.h"
#include "QMCDrivers/tests/SetupPools.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/SampleStack.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Utilities/ProjectData.h"
#include "Utilities/for_testing/Catch2Approx.h"
#include <MinimalHamiltonianPool.h>
#include <MinimalParticlePool.h>
#include <MinimalWaveFunctionPool.h>

namespace qmcplusplus
{
namespace testing
{
class QMCFixedSampleLinearOptimizeInputTest
{
public:
  static bool read(QMCFixedSampleLinearOptimize& driver, xmlNodePtr node) { return driver.m_param.put(node); }

  static bool read(QMCFixedSampleLinearOptimizeBatched& driver, xmlNodePtr node)
  {
    return driver.m_param.put(node);
  }

  static void checkLegacyDefaults(const QMCFixedSampleLinearOptimize& driver)
  {
    CHECK(driver.nblocks == 1);
    CHECK(driver.nolds == 1);
    CHECK(driver.nkept == 1);
    CHECK(driver.num_shifts == 3);
    CHECK(driver.max_param_change == Catch::Approx(0.3));
    CHECK(driver.Max_iterations == 1);
    CHECK(driver.MinMethod == "OneShiftOnly");
  }

  static void checkLegacyLiveValues(const QMCFixedSampleLinearOptimize& driver)
  {
    CHECK(driver.nblocks == 5);
    CHECK(driver.nolds == 4);
    CHECK(driver.nkept == 2);
    CHECK(driver.num_shifts == 6);
    CHECK(driver.max_param_change == Catch::Approx(0.25));
    CHECK(driver.Max_iterations == 7);
    CHECK(driver.MinMethod == "adaptive");
  }

  static void checkBatchedDefaults(const QMCFixedSampleLinearOptimizeBatched& driver)
  {
    CHECK(driver.options_LMY_.nblocks == 1);
    CHECK(driver.options_LMY_.nolds == 1);
    CHECK(driver.options_LMY_.nkept == 1);
    CHECK(driver.options_LMY_.num_shifts == 3);
    CHECK(driver.options_LMY_.max_param_change == Catch::Approx(0.3));
    CHECK_FALSE(driver.options_LMY_.store_samples);
    CHECK(driver.Max_iterations == 1);
    CHECK(driver.MinMethod == "OneShiftOnly");
  }

  static void checkBatchedLiveValues(const QMCFixedSampleLinearOptimizeBatched& driver)
  {
    CHECK(driver.options_LMY_.nblocks == 5);
    CHECK(driver.options_LMY_.nolds == 4);
    CHECK(driver.options_LMY_.nkept == 2);
    CHECK(driver.options_LMY_.num_shifts == 6);
    CHECK(driver.options_LMY_.max_param_change == Catch::Approx(0.25));
    CHECK(driver.options_LMY_.store_samples);
    CHECK(driver.Max_iterations == 7);
    CHECK(driver.MinMethod == "adaptive");
  }
};
} // namespace testing

TEST_CASE("QMCFixedSampleLinearOptimize retired comparison samples", "[drivers][wfopt]")
{
  using testing::QMCFixedSampleLinearOptimizeInputTest;

  Communicate* comm = OHMMS::Controller;
  ProjectData project("test", ProjectData::DriverVersion::LEGACY);
  const SimulationCell simulation_cell;
  MCWalkerConfiguration walkers(simulation_cell);
  walkers.setName("e");
  walkers.create({1});
  walkers.createWalkers(1);
  TrialWaveFunction wavefunction(project.getRuntimeOptions());
  QMCHamiltonian hamiltonian;

  QMCFixedSampleLinearOptimize driver(project, walkers, wavefunction, hamiltonian, comm);
  QMCFixedSampleLinearOptimizeInputTest::checkLegacyDefaults(driver);

  Libxml2Document retired_only_doc;
  REQUIRE(retired_only_doc.parseFromString(R"(<qmc><parameter name="nsamp_comp">23</parameter></qmc>)"));
  CHECK_FALSE(QMCFixedSampleLinearOptimizeInputTest::read(driver, retired_only_doc.getRoot()));
  QMCFixedSampleLinearOptimizeInputTest::checkLegacyDefaults(driver);

  Libxml2Document live_doc;
  REQUIRE(live_doc.parseFromString(R"(
<qmc>
  <parameter name="nsamp_comp">not-an-integer</parameter>
  <parameter name="nblocks">5</parameter>
  <parameter name="nolds">4</parameter>
  <parameter name="nkept">2</parameter>
  <parameter name="num_shifts">6</parameter>
  <parameter name="max_param_change">0.25</parameter>
  <parameter name="max_its">7</parameter>
  <parameter name="MinMethod">adaptive</parameter>
</qmc>)"));
  CHECK(QMCFixedSampleLinearOptimizeInputTest::read(driver, live_doc.getRoot()));
  QMCFixedSampleLinearOptimizeInputTest::checkLegacyLiveValues(driver);
}

TEST_CASE("QMCFixedSampleLinearOptimizeBatched retired comparison samples", "[drivers][wfopt]")
{
  using namespace testing;

  Communicate* comm = OHMMS::Controller;
  ProjectData project("test", ProjectData::DriverVersion::BATCH);
  auto particle_pool = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool =
      MinimalWaveFunctionPool::make_diamondC_1x1x1(project.getRuntimeOptions(), comm, particle_pool);
  auto hamiltonian_pool = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);

  QMCDriverInput qmcdriver_input;
  VMCDriverInput vmcdriver_input;
  WalkerConfigurations walker_confs;
  RandomNumberGeneratorPool rng_pool(1);
  SampleStack samples;
  QMCFixedSampleLinearOptimizeBatched driver(
      project, std::move(qmcdriver_input), std::move(vmcdriver_input), walker_confs,
      MCPopulation(comm->size(), comm->rank(), *particle_pool.getParticleSet("e"),
                   wavefunction_pool.getWaveFunction().value().get(),
                   hamiltonian_pool.getHamiltonian().value().get()),
      rng_pool.getRngRefs(), samples, comm);
  QMCFixedSampleLinearOptimizeInputTest::checkBatchedDefaults(driver);

  Libxml2Document retired_only_doc;
  REQUIRE(retired_only_doc.parseFromString(
      R"(<qmc><parameter name="options_LMY_.nsamp_comp">23</parameter></qmc>)"));
  CHECK_FALSE(QMCFixedSampleLinearOptimizeInputTest::read(driver, retired_only_doc.getRoot()));
  QMCFixedSampleLinearOptimizeInputTest::checkBatchedDefaults(driver);

  Libxml2Document live_doc;
  REQUIRE(live_doc.parseFromString(R"(
<qmc>
  <parameter name="options_LMY_.nsamp_comp">not-an-integer</parameter>
  <parameter name="options_LMY_.nblocks">5</parameter>
  <parameter name="options_LMY_.nolds">4</parameter>
  <parameter name="options_LMY_.nkept">2</parameter>
  <parameter name="options_LMY_.num_shifts">6</parameter>
  <parameter name="options_LMY_.max_param_change">0.25</parameter>
  <parameter name="store_samples">yes</parameter>
  <parameter name="max_its">7</parameter>
  <parameter name="MinMethod">adaptive</parameter>
</qmc>)"));
  CHECK(QMCFixedSampleLinearOptimizeInputTest::read(driver, live_doc.getRoot()));
  QMCFixedSampleLinearOptimizeInputTest::checkBatchedLiveValues(driver);
}
} // namespace qmcplusplus
