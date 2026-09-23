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
#include "Utilities/for_testing/Catch2Approx.h"

#include "Message/Communicate.h"
#include "Message/UniformCommunicateError.h"
#include "QMCDrivers/DMC/DMCDriverInput.h"
#include "QMCDrivers/DMC/DMCBatched.h"
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

TEST_CASE("DMCBatched rejects L2 diffusion with spinors", "[drivers]")
{
  using namespace testing;
  RandomNumberGeneratorPool rng_pool(1);
  ProjectData test_project;
  Communicate* comm = OHMMS::Controller;

  Libxml2Document doc;
  REQUIRE(doc.parseFromString(R"(<qmc method="dmc"><parameter name="L2_diffusion">yes</parameter></qmc>)"));
  xmlNodePtr node = doc.getRoot();
  QMCDriverInput qmcdriver_input;
  qmcdriver_input.readXML(node);
  DMCDriverInput dmcdriver_input;
  dmcdriver_input.readXML(node);

  auto particle_pool = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool =
      MinimalWaveFunctionPool::make_diamondC_1x1x1(test_project.getRuntimeOptions(), comm, particle_pool);
  auto hamiltonian_pool = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);
  particle_pool.getParticleSet("e")->setSpinor(true);

  WalkerConfigurations walker_confs;
  auto construct_driver = [&]() {
    DMCBatched dmcdriver(test_project, std::move(qmcdriver_input), nullptr, std::move(dmcdriver_input), walker_confs,
                         MCPopulation(comm->size(), comm->rank(), *particle_pool.getParticleSet("e"),
                                      wavefunction_pool.getWaveFunction().value(),
                                      hamiltonian_pool.getHamiltonian().value()),
                         rng_pool.getRngRefs(), comm);
  };

  try
  {
    construct_driver();
    FAIL("DMCBatched accepted L2 diffusion with a spinor ParticleSet");
  }
  catch (const UniformCommunicateError& error)
  {
    CHECK(std::string(error.what()) == "L2 diffusion is not supported for spinor particle sets.");
  }
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
  Communicate* comm;
  comm = OHMMS::Controller;
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
                       MCPopulation(comm->size(), comm->rank(), *particle_pool.getParticleSet("e"),
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
  CHECK(dmcdriver.get_num_living_walkers() == 8);
  const QMCTraits::IndexType reserved_walkers = dmcdriver.get_num_living_walkers() + dmcdriver.get_num_dead_walkers();
  CHECK(reserved_walkers == 10);
  // What else should we expect after process
}
#endif

} // namespace qmcplusplus
